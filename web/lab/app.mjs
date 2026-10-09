import {estimateRcs,diameterRoots,snrDb} from './rcs.mjs?v=20261009residual11';
import {initialState,restoreControls,initializeState,saveState} from './gui-state.mjs?v=20261009residual11';
import { Heatmap, LinePlot, db } from "./plots.mjs?v=20261009residual11";
const $ = (id) => document.getElementById(id),
  num = (id) => Number($(id)?.value ?? state[id === "noiseStart" ? "bgStart" : id === "noiseStop" ? "bgStop" : id]),
  worker = new Worker("worker.mjs?v=20261009residual11", { type: "module" });
let meta,
  period,
  framePeriod,
  loaded = false,
  busy = false,
  power,
  viewImages,
  fftPeak = 1,
  result,
  compareResults = [];
let analysisComplete=false;
let restoreFit=initialState?.fit===true;
const state = { start: 0, stop: 6250, chirp: 625, pulses: 8, rx: 0, component:0, bgStart:0,bgStop:125,scanStart:625,scanStop:750 };
function status(text, error = false) {
  $("status").textContent = text;
  $("status").classList.toggle("error", error);
}
function setBusy(value) {
  busy = value;
  window.dispatchEvent(new CustomEvent("fmcw-busy", { detail: value }));
  for (const el of document.querySelectorAll(
    ".settings input,.settings select,.settings button,#matchSettings input,#matchSettings select,#search,#scan",
  ))
    el.disabled = value || !loaded;
  $("cancel").disabled = !value;
  $("download").disabled=value||!analysisComplete;
  for (const id of ["xN", "yN", "vN"]) $(id).disabled = value || !loaded || $("algorithm").value === "fft";
  $("phaseLoss").disabled = value || !loaded || $("algorithm").value !== "fft";
  for (const map of maps) map.drawWindow();
}
function assertRange(a, b, label = "chirp") {
  if (
    !Number.isInteger(a) ||
    !Number.isInteger(b) ||
    a < 0 ||
    b > meta.total_chirps ||
    a >= b
  )
    throw Error(`Enter a valid, nonempty ${label} interval.`);
}
function clock(chirp) {
  const k = Math.round(chirp),
    frame = Math.floor(k / meta.chirps_per_frame),
    within = k % meta.chirps_per_frame;
  const elapsed = frame * framePeriod + within * period + meta.parameters.T_adc;
  return elapsed.toFixed(6);
}
function bands() {
  return [
    {
      start: num("bgStart"),
      stop: num("bgStop"),
      color: "#78b7ff18",
      line: "#78b7ff88",
    },
    {start:num("scanStart"),stop:num("scanStop"),color:"#a855f710",line:"#7e22ce"},
    {
      start: state.chirp,
      stop: state.chirp + state.pulses,
      color: "#f1cf7433",
      line: "#b45c00",
      strong: true,
    },
  ];
}
const select = (k) => {
  if (busy) return;
  selectChirp(Math.floor(k));
};
const maps = [
    new Heatmap("rawPlot", select),
    new Heatmap("subPlot", select),
    new Heatmap("fftPlot", select),
  ],
  vx = new Heatmap("vxPlot"),
  vy = new Heatmap("vyPlot");
const iq = new LinePlot("iqPlot"),
  spectrum = new LinePlot("spectrumPlot"),
  fit = new LinePlot("fitPlot"),
  rangePlot = new LinePlot("rangePlot"),
  residualPlot = new LinePlot("residualPlot"),
  scan = new LinePlot("scanPlot", (k) => selectChirp(Math.round(k)));
const scanPoints=new Map();
const scanPlots=Object.fromEntries(['scanRange','scanVelocity','scanAcceleration','scanPhases','scanSnr','scanRcs','scanDiameter'].map(id=>[id,new LinePlot(id)]));
function updateScan(r) {
  const [r0,v0,a0]=r.best;
  let kinematics=[r0,v0,a0];
  if(r.model==='geometry'){
    const x=r.best[2]*r.midpoint-r.best[0],range=Math.hypot(x,r.best[1]);
    kinematics=[range,r.best[2]*x/range,r.best[2]**2*r.best[1]**2/range**3];
  }
  scanPoints.set(r.start,{start:r.start,midpoint:r.midpoint,kinematics,pulses:r.pulses,best:Array.from(r.best),model:r.model,grid:r.grid,referenceStarts:r.referenceStarts,beam:r.beam});
  drawScanHistory();
}
function drawScanHistory() {
  if(!scanPoints.size)return;
  const points=[...scanPoints.values()].sort((a,b)=>a.start-b.start),x=points.map(p=>Number(clock(p.start))+p.midpoint-meta.parameters.T_adc);
  const cfg={x0:Number(clock(num('scanStart'))),x1:Number(clock(num('scanStop'))),xLabel:'Seconds since file start',mode:'scatter',xFormat:t=>t.toFixed(5)};
  const line=(name,color,y)=>({name,color,x,y});
  const colors=['#0072b2','#d55e00','#009e73','#cc79a7','#111111'];
  for(const [id,index,label] of [['scanRange',0,'r₀ (m)'],['scanVelocity',1,'v₀ (m/s)'],['scanAcceleration',2,'a₀ (10⁶ m/s²)']])scanPlots[id].set([line(label,colors[index],points.map(p=>p.kinematics[index]))],{...cfg,yLabel:label,decimals:2,yFormat:id==='scanAcceleration'?v=>(v/1e6).toFixed(2):undefined});
  const pairs=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]];
  scanPlots.scanPhases.set(pairs.map(([a,b],i)=>line(`RX${b} − RX${a}`,['#0072b2','#d55e00','#009e73','#cc79a7','#111111','#e69f00'][i],points.map(p=>{const d=p.beam.phases[b]-p.beam.phases[a];return Math.atan2(Math.sin(d),Math.cos(d))*180/Math.PI;}))),{...cfg,yLabel:'Phase difference (°)',y0:-180,y1:180});
  scanPlots.scanSnr.set([...Array.from({length:4},(_,i)=>line(`RX${i}`,colors[i],points.map(p=>snrDb(p.beam.single[i])))),line('Four-RX beamformed',colors[4],points.map(p=>snrDb(p.beam.peak)))],{...cfg,mode:'scatter',yLabel:'SNR (dB)'});
  const control=(id,fallback)=>Number(document.getElementById(id)?.value??fallback);
  const frequency=meta.parameters.f_start+meta.parameters.freq_slope*(meta.parameters.T_adc+(meta.samples-1)/(2*meta.parameters.fs));
  const settings={temperature:control('rcsTemperature',9000),frequency,txPowerDbm:control('rcsPower',12),txGainDbi:control('rcsTxGain',6),rxGainDbi:control('rcsRxGain',6),lossDb:control('rcsLoss',0)};
  const minimum=control('rcsDMin',.01)/1000,maximum=control('rcsDMax',20)/1000,key=JSON.stringify([settings,minimum,maximum]);
  const valid=Object.values(settings).every(Number.isFinite)&&settings.temperature>0&&minimum>=1e-6&&maximum>=minimum&&maximum<=.1;
  for(const point of points)if(point.rcsKey!==key){
    point.rcsKey=key;
    point.rcs=valid&&point.beam.noiseCalibrated!==false?[...point.beam.single,point.beam.peak].map((score,i)=>{
      const sigma=estimateRcs(score,{...settings,range:point.kinematics[0],T_coh:point.beam.T_coh??point.beam.effectiveTime},i===4?4:1);
      return {sigma,diameters:diameterRoots(sigma,frequency,minimum,maximum).map(d=>1000*d)};
    }):Array.from({length:5},()=>({sigma:NaN,diameters:[]}));
  }
  scanPlots.scanRcs.set(Array.from({length:5},(_,i)=>line(i===4?'Four-RX beamformed':`RX${i}`,colors[i],points.map(p=>p.rcs[i].sigma>0?10*Math.log10(p.rcs[i].sigma):NaN))),{...cfg,yLabel:'Estimated RCS (dBsm)',decimals:1});
  scanPlots.scanDiameter.set(Array.from({length:5},(_,i)=>{
    const count=Math.max(...points.map(p=>p.rcs[i].diameters.length));
    return Array.from({length:count},(_,root)=>({name:`${i===4?'Four-RX beamformed':`RX${i}`} · solution ${root+1}`,color:colors[i],x,y:points.map(p=>p.rcs[i].diameters[root]??NaN)}));
  }).flat(),{...cfg,mode:'scatter',yLabel:'PEC sphere diameter (mm)',decimals:2});
}
document.addEventListener('change',event=>{if(/^rcs/.test(event.target.id)&&event.target.id!=='rcsRange')drawScanHistory();});
for (const map of maps)
  map.enableWindow({
    getRange: () => [num("bgStart"), num("bgStop")],
    enabled: () => loaded && !busy,
    limit: () => meta.total_chirps,
    preview: (start, stop) => {
      state.bgStart = start;
      state.bgStop = stop;
      maps.forEach((map) => map.overlay({ bands: bands() }));
    },
    commit: () => {
      invalidate();
      status("Updating the complex background mean…");
      applyBackground();
    },
  });
for (const map of maps)
  map.enableWindow({
    label: "Coherent integration",
    getLabel: () => `Coherent integration · ${state.pulses} chirps`,
    moveLabel: "Move coherent integration window",
    tone: "coherent-window",
    resizable: false,
    getRange: () => [state.chirp, state.chirp + state.pulses],
    enabled: () => loaded && !busy,
    limit: () => meta.total_chirps,
    preview: (start) => {
      const end =
        (Math.floor(start / meta.chirps_per_frame) + 1) * meta.chirps_per_frame;
      state.chirp = Math.max(0, Math.min(start, end - state.pulses));
      $("chirp").value = state.chirp;
      $("chirpSlider").value = state.chirp;
      maps.forEach((plot) => plot.overlay({ bands: bands() }));
    },
    commit: () => {
      invalidate();
      selectChirp(state.chirp);
      run();
    },
  });
for(const map of maps)map.enableWindow({
  label:'Analyze interval',tone:'analysis-window',
  getLabel:()=>`Analyze · ${num('scanStart')}–${num('scanStop')-1}`,
  getRange:()=>[num('scanStart'),num('scanStop')],
  enabled:()=>loaded&&!busy,limit:()=>meta.total_chirps,
  preview:(start,stop)=>{state.scanStart=start;state.scanStop=stop;maps.forEach(m=>m.overlay({bands:bands()}));},
  commit:()=>window.dispatchEvent(new Event('fmcw-view-state')),
});
for (const [index,map] of maps.entries()) map.enableZoom({
  enabled:()=>loaded && !busy,
  reset:false,
  onZoom:bounds=>{
    const targets=index<2?maps.slice(0,2):[maps[2]];
    for(const target of targets)target.zoomTo({y0:bounds.y0,y1:bounds.y1});
    state.start=Math.max(0,Math.floor(bounds.x0));
    state.stop=Math.min(meta.total_chirps,Math.ceil(bounds.x1));
    customView();
  },
});
function common() {
  return {
    x0: state.start,
    x1: state.stop,
    y0: meta.parameters.T_adc * 1e6,
    y1: (meta.parameters.T_adc + (meta.samples - 1) / meta.parameters.fs) * 1e6,
    xFormat: (x) => Math.min(state.stop - 1, Math.round(x)).toString(),
    topFormat: (x) => clock(Math.min(state.stop - 1, Math.round(x))),
    topLabel: "Time since file start (s; reconstructed)",
    xLabel: "Concatenated chirp number (frame gaps omitted)",
    yLabel: "Fast time from chirp ramp start (µs)",
    bands: bands(),
    frames: meta.chirps_per_frame,
  };
}
function drawView() {
  saveState();
  if (!viewImages) return;
  const base = common();
  for (let i = 0; i < 3; i++) {
    const image = viewImages[i];
    let values = image.values,
      c = { ...base };
    if (i < 2) {
      const range = num(i === 0 ? "rawScale" : "subScale");
      c = { ...c, lo: -range, hi: range, kind: 0, colorLabel: "ADC counts" };
    } else {
      values = Float32Array.from(values, (x) =>
        Number.isFinite(x) ? db(x / fftPeak) : NaN,
      );
      c = {
        ...c,
        y0: -meta.parameters.fs / 2e6,
        y1: (meta.parameters.fs / 2 - meta.parameters.fs / image.height) / 1e6,
        lo: -num("span"),
        hi: 0,
        kind: 1,
        yLabel: "Stored-I/Q beat frequency (MHz)",
        colorLabel: "dB / raw spectral peak",
      };
    }
    maps[i].set(values, state.stop - state.start, image.height, c);
  }
}
function view() {
  const a = state.start,
    b = state.stop;
  assertRange(a, b);
  state.start = a;
  state.stop = b;
  worker.postMessage({
    type: "view",
    start: a,
    stop: b,
    component: num("component"),
    fftSub: !!num("fftSub"),
  });
}
function customView() {view();}
function invalidate() {
  window.dispatchEvent(new Event("fmcw-invalidated"));
  const hadResult = !!result;
  if($("algorithm").value==='fft') for(const id of ['xN','yN','vN']) $(id).value='';
  result = undefined;
  analysisComplete=false;$("download").disabled=true;
  $("progress").value = 0;
  vx.root.style.display = "none";
  vy.root.style.display = "none";
  fit.root.style.display = "none";
  rangePlot.root.style.display = "none";
  residualPlot.root.style.display = "none";
  scan.root.style.display = "none";
  $("compareResults").innerHTML = "";
  $("matchSummary").textContent =
    hadResult
      ? "Fit cleared; search again."
      : "Search selected train.";
}
function updateSelection() {
  for (const m of maps) m.overlay({ bands: bands() });
  worker.postMessage({ type: "trace", chirp: state.chirp });
}
function selectChirp(k) {
  if (!meta) return;
  k = Math.max(0, Math.min(meta.total_chirps - 1, k));
  const n = num("pulses"),
    end = (Math.floor(k / meta.chirps_per_frame) + 1) * meta.chirps_per_frame;
  if (k + n > end) k = end - n;
  if (k !== state.chirp) invalidate();
  state.pulses = n;
  state.chirp = k;
  $("chirp").value = k;
  $("chirpSlider").value = k;
  updateSelection();
}
function applyTiming() {
  const previous = period;
  const tail = Number(initialState?.controls?.tail ?? meta.ramp_tail_us_assumption) * 1e-6;
  period = Number(initialState?.timing?.period ?? (meta.parameters.T_adc + meta.samples / meta.parameters.fs + meta.parameters.T_idle + tail));
  framePeriod = Number(initialState?.timing?.framePeriod ?? (initialState?.controls?.framePeriod ? Number(initialState.controls.framePeriod)/1000 : meta.chirps_per_frame*period));
  if(!Number.isFinite(period)||period<=0||!Number.isFinite(framePeriod)||framePeriod+1e-12<meta.chirps_per_frame*period)throw Error('Invalid reconstructed timing.');
  if(previous!==undefined&&previous!==period)invalidate();
  worker.postMessage({type:'period',period});
  drawView();
  updateSelection();
}
function trace(msg) {
  if (msg.chirp !== state.chirp) return;
  const n = meta.samples,
    nf = 2 ** Math.ceil(Math.log2(n)),
    a = msg.values,
    x = Float64Array.from(
      { length: n },
      (_, j) => (meta.parameters.T_adc + j / meta.parameters.fs) * 1e6,
    ),
    offset = $("traceSub").checked ? 2 * n : 0;
  const desc = `Global ${msg.chirp} · frame ${Math.floor(msg.chirp / meta.chirps_per_frame)}, chirp ${msg.chirp % meta.chirps_per_frame} · time since file start ${clock(msg.chirp)} s (reconstructed)${msg.valid ? "" : " · FLAGGED PADDING: excluded from filters/spectra"}`;
  $("chirpInfo").textContent = desc;
  iq.set(
    [
      { name: "Re", color: "#b45c00", x, y: a.slice(offset, offset + n) },
      {
        name: "Im",
        color: "#007c78",
        x,
        y: a.slice(offset + n, offset + 2 * n),
      },
    ],
    { x0: x[0], x1: x[n - 1], xLabel: "Fast time (µs)", yLabel: "ADC counts",mode:"line" },
  );
  const fx = Float64Array.from(
    { length: nf },
    (_, j) => ((j - nf / 2) * meta.parameters.fs) / nf / 1e6,
  );
  spectrum.set(
    [
      {
        name: "Original",
        color: "#555555",
        x: fx,
        y: Float32Array.from(a.slice(4 * n, 4 * n + nf), (v) =>
          db(v / fftPeak),
        ),
      },
      {
        name: "Quiet mean removed",
        color: "#007c78",
        x: fx,
        y: Float32Array.from(a.slice(4 * n + nf), (v) => db(v / fftPeak)),
      },
    ],
    {
      x0: fx[0],
      x1: fx[nf - 1],
      y0: -num("span"),
      y1: 5,
      xLabel: "Stored-I/Q beat frequency (MHz)",
      mode:"line",
      yLabel: "dB / raw spectral peak",
    },
  );
}
function grid() {
  const g = Object.fromEntries(
    ["xMin", "xMax", "xN", "yMin", "yMax", "yN", "vMin", "vMax", "vN"].map(
      (k) => [k, num(k)],
    ),
  );
  for (const a of ["x", "y", "v"])
    if (
      !Number.isFinite(g[a + "Min"]) ||
      !Number.isFinite(g[a + "Max"]) ||
      g[a + "Min"] > g[a + "Max"] ||
      ($("algorithm").value === "direct" && (!Number.isInteger(g[a + "N"]) || g[a + "N"] < 1))
    )
      throw Error("Invalid trajectory grid.");
  if (($("algorithm").value === "fft" ? g.xMin <= 0 : g.yMin <= 0) || Math.max(Math.abs(g.vMin), Math.abs(g.vMax)) > 100000)
    throw Error("Use positive range (radial FFT) or y₀ (direct geometry), and |v| ≤ 100 km/s.");
  const cells = g.xN * g.yN * g.vN;
  if ($("algorithm").value === "direct" && cells > 150000)
    throw Error(
      "Limit the bank to 150,000 templates; use a coarse search then refine.",
    );
  return g;
}
function run(compare = false, timeScan = false) {
  if (busy) return;
  try {
    applyTiming();
    assertRange(num("bgStart"), num("bgStop"), "quiet-mean chirp");
    const g = grid();
    const loss = num("phaseLoss");
    if ($("algorithm").value === "fft" && (!Number.isFinite(loss) || loss <= 0 || loss > 50))
      throw Error("Choose a correction-phase loss greater than 0% and at most 50%.");
    assertRange(num("noiseStart"), num("noiseStop"), "noise-reference chirp");
    const needed = compare ? 16 : num("pulses");
    if (
      state.chirp + needed >
      (Math.floor(state.chirp / meta.chirps_per_frame) + 1) *
        meta.chirps_per_frame
    )
      if(!timeScan)throw Error("Select an earlier chirp in this frame so all requested pulses fit.");
    const job = {
      type: "search",
      start: state.chirp,
      pulses: num("pulses"),
      period,
      receiver: $("receiver").checked,
      noiseStart: num("noiseStart"),
      noiseStop: num("noiseStop"),
      grid: g,
      algorithm: $("algorithm").value,
      computeBackend: $("computeBackend").value,
      loss,
      compare,
      scan: timeScan,
      scanStart: num("scanStart"),
      scanStop: num("scanStop"),
      scanStride: 1,
      beamSteps: Number(document.getElementById("beamSteps")?.value??10),
    };
    if (timeScan) {
      assertRange(job.scanStart, job.scanStop, "time-scan chirp");
      if(job.scanStop-job.scanStart<job.pulses)throw Error("Analyze interval must contain a complete coherent train.");
      if (!Number.isInteger(job.scanStride) || job.scanStride < 1)
        throw Error("Start stride must be a positive integer.");

    }
    analysisComplete=false;
    setBusy(true);
    compareResults = [];
    $("compareResults").innerHTML = "";
    $("progress").value = 0;
    status(job.algorithm === "fft" ? "Building the explicit r₀ / v₀ / a₀ grid…" : `Searching ${g.xN * g.yN * g.vN} trajectory templates…`);
    // Queue the background only after validation: its asynchronous reply must
    // not replace a rejected search's error with a successful-update message.
    worker.postMessage({
      type: "background",
      start: num("bgStart"),
      stop: num("bgStop"),
    });
    worker.postMessage(job);
  } catch (e) {
    status(e.message, true);
  }
}
function axis(lo, hi, n) {
  return Float64Array.from({ length: n }, (_, i) =>
    n === 1 ? lo : lo + ((hi - lo) * i) / (n - 1),
  );
}
function showMatch(r) {
  window.dispatchEvent(new CustomEvent("fmcw-match", { detail: { meta, match: r } }));
  result = r;
  if(r.scanPoint)updateScan(r);
  vx.resetZoom();vy.resetZoom();
  state.pulses = r.pulses;
  state.chirp = r.start;
  $("chirp").value = r.start;
  $("chirpSlider").value = r.start;
  if(r.scanPoint)updateSelection();
  vx.root.style.display = "block";
  vy.root.style.display = "block";
  maps.forEach((m) => m.overlay({ bands: bands() }));
  saveState();
  const g = r.grid,
    nv = r.radialSpec ? r.radialSpec[2] : r.cells ? Math.max(1,Math.min(48,Math.ceil(Math.sqrt(r.cube.length)))) : g.vN,
    nx = r.radialSpec ? r.radialSpec[0] : r.cells ? nv : g.xN,
    ny = r.radialSpec ? r.radialSpec[1] : r.cells ? nv : g.yN,
    a = new Float32Array(nx * nv).fill(-Infinity),
    b = new Float32Array(ny * nv).fill(-Infinity);
  if (r.cells) {
    const project = (dest, lo, hi, axis, rows) => {
      for (let cell=0;cell<r.cube.length;cell++) {
        const score=snrDb(r.cube[cell]), point=r.nodes.subarray(cell*3,cell*3+3);
        const row=Math.min(rows-1,Math.max(0,Math.floor((point[axis]-lo)/(hi-lo || 1)*rows)));
        const col=Math.min(nv-1,Math.max(0,Math.floor((point[2]-g.vMin)/(g.vMax-g.vMin || 1)*nv)));
        dest[row*nv+col]=Math.max(dest[row*nv+col],score);
      }
    };
    project(a,g.xMin,g.xMax,0,nx); project(b,g.yMin,g.yMax,1,ny);
  } else for (let ix = 0; ix < nx; ix++)
    for (let iy = 0; iy < ny; iy++)
      for (let iv = 0; iv < nv; iv++) {
        const score = snrDb(r.cube[(ix * ny + iy) * nv + iv]);
        if (score > a[ix * nv + iv]) a[ix * nv + iv] = score;
        if (score > b[iy * nv + iv]) b[iy * nv + iv] = score;
      }
  let gridPeakIndex=0, gridPeak=-Infinity;
  for(let i=0;i<r.cube.length;i++) if(r.cube[i]>gridPeak){gridPeak=r.cube[i];gridPeakIndex=i;}
  const markerPoint=r.radialSpec ? [g.xMin+(g.xMax-g.xMin)*Math.floor(gridPeakIndex/(ny*nv))/(nx-1 || 1),g.yMin+(g.yMax-g.yMin)*(Math.floor(gridPeakIndex/nv)%ny)/(ny-1 || 1),g.vMin+(g.vMax-g.vMin)*(gridPeakIndex%nv)/(nv-1 || 1)] : [r.best[0],r.best[1],r.best[2]];
  const peak = snrDb(r.radialSpec ? gridPeak : r.best[3]),
    low = Math.min(0, peak - 35),
    high = peak;
  const opts = {
    x0: g.vMin,
    x1: g.vMax === g.vMin ? g.vMax + 1 : g.vMax,
    lo: low,
    hi: high,
    topLabel: `Peak matched-filter value = ${peak.toFixed(2)} dB`,
    colorDecimals: 2,
    kind: 1,
    xLabel: r.radialSpec ? "Radial velocity v₀ (m/s; positive receding)" : "Along-track velocity v (m/s)",
    xMath:r.radialSpec?String.raw`v_0\;(\mathrm{m\,s^{-1}})`:String.raw`v\;(\mathrm{m\,s^{-1}})`,
    xFormat: (x) => x.toFixed(0),
    colorLabel: "SNR (dB)",
  };
  vx.set(a, nv, nx, {
    ...opts,
    y0: g.xMin,
    y1: g.xMax === g.xMin ? g.xMax + 0.01 : g.xMax,
    yLabel: r.radialSpec ? "Midpoint range r₀ (m)" : "Along-track offset x₀ (m)",
    yMath:r.radialSpec?String.raw`r_0\;(\mathrm{m})`:String.raw`x_0\;(\mathrm{m})`,
    yFormat: (x) => x.toFixed(3),
    marker: {
      column: Math.round(
        ((markerPoint[2] - g.vMin) / (g.vMax - g.vMin || 1)) * (nv - 1),
      ),
      row: Math.round(
        ((markerPoint[0] - g.xMin) / (g.xMax - g.xMin || 1)) * (nx - 1),
      ),
    },
  });
  vy.set(b, nv, ny, {
    ...opts,
    y0: g.yMin,
    y1: g.yMax === g.yMin ? g.yMax + 0.01 : g.yMax,
    yLabel: r.radialSpec ? "a₀ (10⁶ m/s²)" : "y₀ (m)",
    yMath:r.radialSpec?String.raw`a_0\;(10^6\,\mathrm{m\,s^{-2}})`:String.raw`y_0\;(\mathrm{m})`,
    yFormat: (x) => r.radialSpec ? (x/1e6).toFixed(2) : x.toFixed(2),
    marker: {
      column: Math.round(
        ((markerPoint[2] - g.vMin) / (g.vMax - g.vMin || 1)) * (nv - 1),
      ),
      row: Math.round(
        ((markerPoint[1] - g.yMin) / (g.yMax - g.yMin || 1)) * (ny - 1),
      ),
    },
  });
  const [x, y, v, , , noiseCount] = r.best;
  const overlap =
    r.start < num("bgStop") && r.start + r.pulses > num("bgStart");
  if(r.radialSpec) { $("xN").value=nx; $("yN").value=ny; $("vN").value=nv; }
  const boundsHit=r.radialSpec&&[[r.best[0],g.xMin,g.xMax],[r.best[1],g.vMin,g.vMax],[r.best[2],g.yMin,g.yMax]].some(([v,lo,hi])=>hi>lo&&Math.min(Math.abs(v-lo),Math.abs(v-hi))<1e-5*(hi-lo));
  $("matchSummary").textContent =
    `${r.searchReceivers??1} RX incoherent · ${r.pulses} chirps (${r.start}–${r.start+r.pulses-1}) · ${snrDb(r.best[3]).toFixed(2)} dB · ${r.radialSpec?`r₀ ${x.toFixed(4)} m, v₀ ${y.toFixed(2)} m/s, a₀ ${v.toPrecision(5)} m/s²`:`x₀ ${x.toFixed(4)} m, y₀ ${y.toFixed(4)} m, v ${v.toFixed(2)} m/s`} · ${r.backend??'CPU'} · ${r.seconds.toFixed(2)} s${r.fallbackReason?' · CPU fallback':''}${overlap?' · overlaps background':''}${boundsHit?' · fit reaches search bound':''}`;
  const xx = [],
    obs = [],
    model = [],
    residualImag = [];
  for (let k = 0; k < r.pulses; k++) {
    for (let j = 0; j < meta.samples; j++) {
      xx.push(
        (k * r.period + meta.parameters.T_adc + j / meta.parameters.fs) * 1e6,
      );
      obs.push(r.observed[2 * (k * meta.samples + j)]);
      model.push(r.fit[2 * (k * meta.samples + j)]);
      residualImag.push(
        r.observed[2 * (k * meta.samples + j) + 1] -
          r.fit[2 * (k * meta.samples + j) + 1],
      );
    }
    if (k < r.pulses - 1) {
      xx.push(
        (k * r.period +
          meta.parameters.T_adc +
          meta.samples / meta.parameters.fs) *
          1e6,
      );
      obs.push(NaN);
      model.push(NaN);
      residualImag.push(NaN);
    }
  }
  fit.set(
    [
      { name: "Re: measured residual", color: "#555555", x: xx, y: obs },
      {
        name: "Re: fitted template",
        color: "#b45c00",
        x: xx,
        y: model,
        width: 1.7,
      },
    ],
    {
      x0: xx[0],
      x1: xx.at(-1),
      xLabel: "Time from train ramp start (µs; gaps are unsampled)",
      yLabel: "ADC counts",mode:"line",
    },
  );
  residualPlot.set(
    [
      {
        name: "Residual Re: data − fit",
        color: "#b45c00",
        x: xx,
        y: Float64Array.from(obs, (value, i) => value - model[i]),
      },
      { name: "Residual Im", color: "#007c78", x: xx, y: residualImag },
    ],
    {
      x0: xx[0],
      x1: xx.at(-1),
      xLabel: "Time from train ramp start (µs; gaps are unsampled)",
      yLabel: "Residual (ADC counts)",mode:"line",
    },
  );
  const rangeTime = Float64Array.from(
    { length: 301 },
    (_, i) => xx[0] + ((xx.at(-1) - xx[0]) * i) / 300,
  );
  const ranges = Float64Array.from(rangeTime, (t) =>
    r.radialSpec ? x+y*(t*1e-6-r.midpoint)+0.5*v*(t*1e-6-r.midpoint)**2 : Math.hypot(v * t * 1e-6 - x, y),
  );
  const minRange = Math.min(...ranges),
    maxRange = Math.max(...ranges);
  const rangePad = Math.max(0.0001, (maxRange - minRange) * 0.08);
  rangePlot.set(
    [
      {
        name: "Best-fit slant range R(t)",
        color: "#b45c00",
        x: rangeTime,
        y: ranges,
        width: 1.7,
      },
    ],
    {
      x0: xx[0],
      x1: xx.at(-1),
      y0: minRange - rangePad,
      y1: maxRange + rangePad,
      xLabel: "Time from train ramp start (µs; continuous trajectory)",
      yLabel: "Slant range (m)",mode:"line",
      decimals: 3,
    },
  );
  if (r.scanResults) {
    const s = r.scanResults,
      starts = [],
      ys = [];
    for (let i = 0; i < s.length; i += 5) {
      starts.push(s[i]);
      ys.push(snrDb(s[i + 1]));
    }
    scan.set(
      [
        {
          name: r.radialSpec ? "Time scan: verified radial FFT candidates" : "Time scan: maximum over x₀, y₀, v and I/Q orientation",
          color: "#007c78",
          x: starts,
          y: ys,
        },
      ],
      {
        x0: starts[0],
        x1: starts.at(-1) + 1,
        xLabel: "Train start chirp (click to select)",
        yLabel: "SNR (dB)",
        decimals: 1,
      },
    );
    state.chirp = r.start;
    $("chirp").value = r.start;
    $("chirpSlider").value = r.start;
    updateSelection();
  }
}
function comparison(results) {
  const radial=!!results[0].radialSpec;
  $("compareResults").innerHTML =
    `<table><thead><tr><th>Pulses</th><th>Sampled µs</th><th>Verified peak dB</th><th>${radial ? "r₀ (m)" : "x₀ (m)"}</th><th>${radial ? "v₀ (m/s)" : "y₀ (m)"}</th><th>${radial ? "a₀ (m/s²)" : "v (m/s)"}</th></tr></thead><tbody>` +
    results
      .map(
        (r) =>
          `<tr><td>${r.pulses}</td><td>${(((r.pulses * meta.samples) / meta.parameters.fs) * 1e6).toFixed(2)}</td><td>${snrDb(r.best[3]).toFixed(2)}</td><td>${r.best[0].toFixed(4)}</td><td>${r.best[1].toFixed(4)}</td><td>${r.best[2].toFixed(2)}</td></tr>`,
      )
      .join("") +
    "</tbody></table>";
  const best = results.reduce((a, b) => (a.best[3] > b.best[3] ? a : b));
  showMatch(best);
  $("pulses").value = best.pulses;
  maps.forEach((m) => m.overlay({ bands: bands() }));
  saveState();
}
worker.onmessage = ({ data: m }) => {
  if (m.type === "loaded") {
    loaded = true;
    status(`RX${m.rx}: ${meta.total_chirps} chirps loaded; SHA-256 verified.`);
    setBusy(false);
  } else if (m.type === "background") {
    invalidate();
    power = m.power;
    if (!busy) status("Complex background mean updated; plots refreshed.");
    applyTiming();
    view();
    updateSelection();
    if(restoreFit){restoreFit=false;setTimeout(()=>run(),0);}
  } else if (m.type === "view") {
    if (m.start !== state.start || m.stop !== state.stop) return;
    viewImages = m.images;
    fftPeak = m.fft_peak;
    drawView();
    updateSelection();
  } else if (m.type === "trace") trace(m);
  else if (m.type === "progress") {
    $("progress").value = m.fraction;
    status(
      `${m.phase}: ${m.pulses} pulses · ${(m.fraction * 100).toFixed(0)}%${m.windows ? " · " + m.windows + " time starts" : ""}`,
    );
  } else if(m.type==='scan-start') {
    scanPoints.clear();$('scanTimeline').hidden=false;
    for(const plot of Object.values(scanPlots))plot.set([],{x0:Number(clock(m.start)),x1:Number(clock(m.stop)),xLabel:'Seconds since file start'});
  } else if(m.type==='scan-progress') {
    $('progress').value=m.fraction;status(`Analyzed train ${m.start} · ${(100*m.fraction).toFixed(1)}%`);
  } else if (m.type === "match") {
    showMatch(m.result);
  } else if (m.type === "complete") {
    analysisComplete=true;
    setBusy(false);
    compareResults = m.results;
    if (m.results.length > 1) comparison(m.results);
    $("progress").value = 1;
    status(
      "Search complete.",
    );
  } else if (m.type === "cancelled") {
    setBusy(false);
    status("Search cancelled; completed plots retained.");
  } else if (m.type === "error") {
    setBusy(false);
    status(m.message, true);
    console.error(m.stack);
  }
};
function safe(fn) {
  return () => {
    try {
      fn();
    } catch (e) {
      status(e.message, true);
    }
  };
}
function applyBackground() {
  assertRange(num('bgStart'),num('bgStop'));
  worker.postMessage({type:'background',start:num('bgStart'),stop:num('bgStop')});
}
$("component").onchange=safe(view);
$("rx").onchange=()=>{
  invalidate();setBusy(true);loaded=false;state.rx=num('rx');
  worker.postMessage({type:'rx',rx:state.rx,bgStart:state.bgStart,bgStop:state.bgStop});
};
$("fftSub").onchange = safe(view);
for (const id of ["rawScale", "subScale", "span"])
  $(id).onchange = safe(drawView);
$("chirp").onchange = () => selectChirp(num("chirp"));
$("chirpSlider").oninput = () => selectChirp(num("chirpSlider"));
$("pulses").onchange = () => {
  invalidate();
  selectChirp(state.chirp);
};
$("traceSub").onchange = updateSelection;
$("full").onclick = safe(() => {
  for(const plot of [...maps,vx,vy])plot.resetZoom();
  state.start=0;state.stop=meta.total_chirps;view();
});
$("scan").onclick = () => run(false, true);
$("cancel").onclick = () => worker.postMessage({ type: "cancel" });
for (const id of [
  "xMin",
  "xMax",
  "xN",
  "yMin",
  "yMax",
  "yN",
  "vMin",
  "vMax",
  "vN",
  "receiver",
  "phaseLoss",
  "computeBackend",
])
  $(id).addEventListener("change", invalidate);
const savedBounds={fft:[.001,3,0,1000000,0,900],direct:[-.6,.6,.05,.3,200,500]};
let activeMethod='fft';
function methodUI() {
 const radial=$("algorithm").value==='fft';
 ['xN','yN','vN'].forEach((id,i)=>$(id).value=radial?'':[31,11,41][i]);
 $("axisR").textContent=radial?'r₀ (m)':'x₀ (m)';$("axisA").textContent=radial?'a₀ (m/s²)':'y₀ (m)';$("axisV").textContent=radial?'v₀ (m/s)':'v (m/s)';
 $("yMin").min=radial?'':'.001';$("yMax").min=radial?'':'.001';$("yMin").step=$("yMax").step=radial?'1000':'.01';
 $("mapRTitle").textContent=radial?'v₀ × r₀ · maximum over a₀':'v × x₀ · maximum over y₀';$("mapATitle").textContent=radial?'v₀ × a₀ · maximum over r₀':'v × y₀ · maximum over x₀';
 $("pulses").options[0].disabled=radial;if(radial && num('pulses')===1) $("pulses").value=8;
 $("modelNote").textContent=radial?'r(t) = r₀ + v₀ h + ½a₀h²; h is seconds from train midpoint. Positive v₀ = receding.':'R(t) = √((v t − x₀)² + y₀²); x₀ along track, y₀ perpendicular.';
 $("gridNote").textContent=radial?'Automatic grid; full bounds retained.':'Uniform grid over the specified bounds.';
 $("mapNote").textContent='Four-RX incoherent matched energy. MAX projections; cross = grid peak.';
}
$("algorithm").onchange = () => {
 savedBounds[activeMethod]=['xMin','xMax','yMin','yMax','vMin','vMax'].map(num);activeMethod=$("algorithm").value;
 ['xMin','xMax','yMin','yMax','vMin','vMax'].forEach((id,i)=>$(id).value=savedBounds[activeMethod][i]);
 invalidate();methodUI();setBusy(false);
};
restoreControls();methodUI();restoreControls();activeMethod=$("algorithm").value;
try {
  setBusy(false);
  const response = await fetch("datasets/test63/metadata.json", {
    cache: "no-store",
  });
  if (!response.ok) throw Error("Recording metadata unavailable");
  meta = await response.json();
  $("record").textContent =
    `Test 63 · ${meta.diameter_mm} mm ${meta.material} ball · ${meta.parameters.speed.toFixed(2)} m/s · ${(meta.parameters.f_start / 1e9).toFixed(1)} GHz · ${meta.frames} frames × ${meta.chirps_per_frame} chirps × ${meta.samples} samples × ${meta.receivers} RX`;
  $("chirpSlider").max = meta.total_chirps - 1;
  $("fftNote").textContent =
    `Complex Hann-window FFT: ${meta.samples} acquired samples, ${2 ** Math.ceil(Math.log2(meta.samples))} FFT points. Fourier resolution ${(meta.parameters.fs / meta.samples / 1000).toFixed(2)} kHz; zero-padded bin spacing ${(meta.parameters.fs / 2 ** Math.ceil(Math.log2(meta.samples)) / 1000).toFixed(2)} kHz. Frequency = range beat + Doppler.`;
  restoreControls();
  state.rx=Number(initialState?.controls?.rx??0);state.chirp=num('chirp');state.pulses=num('pulses');
  state.start=Number(initialState?.view?.start??initialState?.controls?.viewStart??0);
  state.stop=Number(initialState?.view?.stop??(initialState?.controls?.viewStop!==undefined?Number(initialState.controls.viewStop)+1:meta.total_chirps));
  for(const [prefix,a,b] of [['background','bgStart','bgStop'],['analysis','scanStart','scanStop']]){
    state[a]=Number(initialState?.[prefix]?.start??initialState?.controls?.[a]??state[a]);
    state[b]=Number(initialState?.[prefix]?.stop??initialState?.controls?.[b]??state[b]);
  }
  initializeState(()=>({view:{start:state.start,stop:state.stop},background:{start:state.bgStart,stop:state.bgStop},analysis:{start:state.scanStart,stop:state.scanStop},timing:{period,framePeriod},fit:!!result}));
  worker.postMessage({
    type: "init",
    rx:state.rx,
    meta,
    bgStart: num("bgStart"),
    bgStop: num("bgStop"),
  });
} catch (e) {
  status(e.message, true);
  console.error(e);
}
$('download').onclick=async()=>{
  if(busy||!analysisComplete)return;
  const button=$('download');button.disabled=true;
  try{
    drawScanHistory();
    const points=[...scanPoints.values()].sort((a,b)=>a.start-b.start).map(p=>({...p,time:Number(clock(p.start))+p.midpoint-meta.parameters.T_adc}));
    const settings={controls:Object.fromEntries([...document.querySelectorAll('input[id],select[id]')].map(el=>[el.id,el.type==='checkbox'?el.checked:el.value])),background:[state.bgStart,state.bgStop],analysis:[state.scanStart,state.scanStop],period,framePeriod};
    const exporter=new Worker('export-worker.mjs?v=20261009residual11',{type:'module'});
    const bytes=await new Promise((resolve,reject)=>{
      exporter.onmessage=({data})=>{exporter.terminate();data.error?reject(Error(data.error)):resolve(data.bytes);};
      exporter.onerror=e=>{exporter.terminate();reject(Error(e.message));};
      exporter.postMessage({meta,points,result,settings,url:location.href});
    });
    const url=URL.createObjectURL(new Blob([bytes],{type:'application/x-hdf5'})),link=document.createElement('a');
    link.href=url;link.download='fmcw-test63-analysis.h5';link.click();setTimeout(()=>URL.revokeObjectURL(url),30000);
  }catch(error){status(`HDF5 export: ${error.message}`,true);}
  finally{button.disabled=busy||!analysisComplete;}
};
