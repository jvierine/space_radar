import { Heatmap, LinePlot, db } from "./plots.mjs?v=20261009brush2";
const $ = (id) => document.getElementById(id),
  num = (id) => Number($(id).value),
  worker = new Worker("worker.mjs?v=20261009brush2", { type: "module" });
let meta,
  period,
  framePeriod,
  loaded = false,
  busy = false,
  power,
  viewImages,
  fftPeak = 1,
  result,
  compareResults = [],
  clockBase,
  zoomTimer,
  centeredZoom = true;
const state = { start: 0, stop: 6250, chirp: 625, pulses: 8, rx: 0 };
function status(text, error = false) {
  $("status").textContent = text;
  $("status").classList.toggle("error", error);
}
function setBusy(value) {
  busy = value;
  window.dispatchEvent(new CustomEvent("fmcw-busy", { detail: value }));
  for (const el of document.querySelectorAll(
    ".settings input,.settings select,.settings button,#matchSettings input,#matchSettings select,#search,#compare,#scan",
  ))
    el.disabled = value || !loaded;
  $("cancel").disabled = !value;
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
  return (clockBase / 1000 + elapsed).toFixed(6);
}
function bands() {
  return [
    {
      start: num("bgStart"),
      stop: num("bgStop"),
      color: "#78b7ff18",
      line: "#78b7ff88",
    },
    {
      start: num("noiseStart"),
      stop: num("noiseStop"),
      color: "#57d8c013",
      line: "#57d8c088",
    },
    {
      start: state.chirp,
      stop: state.chirp + state.pulses,
      color: "#f1cf7433",
      line: "#f1cf74",
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
for (const map of maps)
  map.enableWindow({
    getRange: () => [num("bgStart"), num("bgStop")],
    enabled: () => loaded && !busy,
    limit: () => meta.total_chirps,
    preview: (start, stop) => {
      $("bgStart").value = start;
      $("bgStop").value = stop;
      maps.forEach((map) => map.overlay({ bands: bands() }));
    },
    commit: () => {
      invalidate();
      status("Updating the complex background mean…");
      $("background").click();
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
      const repeatSearch = !!result;
      invalidate();
      selectChirp(state.chirp);
      if (repeatSearch) run();
    },
  });
for (const [index,map] of maps.entries()) map.enableZoom({
  enabled:()=>loaded && !busy,
  reset:false,
  onZoom:bounds=>{
    const targets=index<2?maps.slice(0,2):[maps[2]];
    for(const target of targets)target.zoomTo({y0:bounds.y0,y1:bounds.y1});
    $('viewStart').value=Math.max(0,Math.floor(bounds.x0));
    $('viewStop').value=Math.min(meta.total_chirps-1,Math.ceil(bounds.x1)-1);
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
    topLabel: "Unix epoch seconds (reconstructed; UTC assumed)",
    xLabel: "Concatenated chirp number (frame gaps omitted)",
    yLabel: "Fast time from chirp ramp start (µs)",
    bands: bands(),
    frames: meta.chirps_per_frame,
  };
}
function drawView() {
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
  const a = num("viewStart"),
    b = num("viewStop") + 1;
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
function applyCenteredZoom(allowBusy = false) {
  if (!meta || (busy && !allowBusy)) return;
  centeredZoom = true;
  const minimum = Math.max(1, state.pulses);
  const factor = Math.pow(meta.total_chirps / minimum, num("zoom") / 100);
  const width = Math.max(minimum, Math.round(meta.total_chirps / factor));
  const center = state.chirp + state.pulses / 2;
  const start = Math.max(
    0,
    Math.min(meta.total_chirps - width, Math.round(center - width / 2)),
  );
  $("viewStart").value = start;
  $("viewStop").value = start + width - 1;
  $("zoomLabel").textContent =
    width === meta.total_chirps
      ? `Full recording · ${width} chirps`
      : `${(meta.total_chirps / width).toFixed(1)}× · ${width} chirps`;
  view();
}
function customView() {
  clearTimeout(zoomTimer);
  centeredZoom = false;
  $("zoom").value = 0;
  view();
  $("zoomLabel").textContent =
    `Custom view · ${state.stop - state.start} chirps`;
}
function invalidate() {
  window.dispatchEvent(new Event("fmcw-invalidated"));
  const hadResult = !!result;
  if($("algorithm").value==='fft') for(const id of ['xN','yN','vN']) $(id).value='';
  result = undefined;
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
      ? "Previous fit cleared because the window or filter settings changed. Click Search selected train to calculate the new fit."
      : "No fit calculated for this train. Click Search selected train.";
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
  $("scanStart").value = k;
  $("scanStop").value = end;
  updateSelection();
  if (centeredZoom && num("zoom") > 0) applyCenteredZoom();
}
function applyTiming() {
  const previous = period;
  const tail = num("tail") * 1e-6;
  period =
    meta.parameters.T_adc +
    meta.samples / meta.parameters.fs +
    meta.parameters.T_idle +
    tail;
  framePeriod = num("framePeriod") / 1000;
  if (
    !Number.isFinite(tail) ||
    tail < 0 ||
    !Number.isFinite(framePeriod) ||
    framePeriod + 1e-12 < meta.chirps_per_frame * period
  )
    throw Error(
      "Frame interval must be at least the complete chirp-train duration.",
    );
  if (previous !== undefined && previous !== period) invalidate();
  worker.postMessage({ type: "period", period });
  $("timingNote").textContent =
    `${(period * 1e6).toFixed(3)} µs assumed chirp period · ${(framePeriod * 1e3).toFixed(5)} ms configured frame interval. Only the recording start clock is stored (timezone unspecified); individual chirp clocks are reconstructed. Ramp tail and frame gaps are not measured in this file.`;
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
  const desc = `Global ${msg.chirp} · frame ${Math.floor(msg.chirp / meta.chirps_per_frame)}, chirp ${msg.chirp % meta.chirps_per_frame} · epoch seconds ${clock(msg.chirp)} (UTC assumed)${msg.valid ? "" : " · FLAGGED PADDING: excluded from filters/spectra"}`;
  $("chirpInfo").textContent = desc;
  iq.set(
    [
      { name: "Re", color: "#f1cf74", x, y: a.slice(offset, offset + n) },
      {
        name: "Im",
        color: "#57d8c0",
        x,
        y: a.slice(offset + n, offset + 2 * n),
      },
    ],
    { x0: x[0], x1: x[n - 1], xLabel: "Fast time (µs)", yLabel: "ADC counts" },
  );
  const fx = Float64Array.from(
    { length: nf },
    (_, j) => ((j - nf / 2) * meta.parameters.fs) / nf / 1e6,
  );
  spectrum.set(
    [
      {
        name: "Original",
        color: "#a2b6c2",
        x: fx,
        y: Float32Array.from(a.slice(4 * n, 4 * n + nf), (v) =>
          db(v / fftPeak),
        ),
      },
      {
        name: "Quiet mean removed",
        color: "#57d8c0",
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
      throw Error(
        "Select an earlier chirp in this frame so all requested pulses fit.",
      );
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
      loss,
      compare,
      scan: timeScan,
      scanStart: num("scanStart"),
      scanStop: num("scanStop"),
      scanStride: num("scanStride"),
    };
    if (timeScan) {
      assertRange(job.scanStart, job.scanStop, "time-scan chirp");
      if (!Number.isInteger(job.scanStride) || job.scanStride < 1)
        throw Error("Start stride must be a positive integer.");
      if ((job.scanStop - job.scanStart) / job.scanStride > 512)
        throw Error("Use at most 512 time starts per scan; increase stride.");
    }
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
  vx.resetZoom();vy.resetZoom();
  state.pulses = r.pulses;
  state.chirp = r.start;
  $("chirp").value = r.start;
  $("chirpSlider").value = r.start;
  vx.root.style.display = "block";
  vy.root.style.display = "block";
  maps.forEach((m) => m.overlay({ bands: bands() }));
  const g = r.grid,
    nv = r.radialSpec ? r.radialSpec[2] : r.cells ? Math.max(1,Math.min(48,Math.ceil(Math.sqrt(r.cube.length)))) : g.vN,
    nx = r.radialSpec ? r.radialSpec[0] : r.cells ? nv : g.xN,
    ny = r.radialSpec ? r.radialSpec[1] : r.cells ? nv : g.yN,
    a = new Float32Array(nx * nv).fill(-Infinity),
    b = new Float32Array(ny * nv).fill(-Infinity);
  if (r.cells) {
    const project = (dest, lo, hi, axis, rows) => {
      for (let cell=0;cell<r.cube.length;cell++) {
        const score=db(r.cube[cell]), point=r.nodes.subarray(cell*3,cell*3+3);
        const row=Math.min(rows-1,Math.max(0,Math.floor((point[axis]-lo)/(hi-lo || 1)*rows)));
        const col=Math.min(nv-1,Math.max(0,Math.floor((point[2]-g.vMin)/(g.vMax-g.vMin || 1)*nv)));
        dest[row*nv+col]=Math.max(dest[row*nv+col],score);
      }
    };
    project(a,g.xMin,g.xMax,0,nx); project(b,g.yMin,g.yMax,1,ny);
  } else for (let ix = 0; ix < nx; ix++)
    for (let iy = 0; iy < ny; iy++)
      for (let iv = 0; iv < nv; iv++) {
        const score = db(r.cube[(ix * ny + iy) * nv + iv]);
        if (score > a[ix * nv + iv]) a[ix * nv + iv] = score;
        if (score > b[iy * nv + iv]) b[iy * nv + iv] = score;
      }
  let gridPeakIndex=0, gridPeak=-Infinity;
  for(let i=0;i<r.cube.length;i++) if(r.cube[i]>gridPeak){gridPeak=r.cube[i];gridPeakIndex=i;}
  const markerPoint=r.radialSpec ? [g.xMin+(g.xMax-g.xMin)*Math.floor(gridPeakIndex/(ny*nv))/(nx-1 || 1),g.yMin+(g.yMax-g.yMin)*(Math.floor(gridPeakIndex/nv)%ny)/(ny-1 || 1),g.vMin+(g.vMax-g.vMin)*(gridPeakIndex%nv)/(nv-1 || 1)] : [r.best[0],r.best[1],r.best[2]];
  const peak = db(r.radialSpec ? gridPeak : r.best[3]),
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
    xFormat: (x) => x.toFixed(0),
    colorLabel: "Matched energy / quiet noise (dB)",
  };
  vx.set(a, nv, nx, {
    ...opts,
    y0: g.xMin,
    y1: g.xMax === g.xMin ? g.xMax + 0.01 : g.xMax,
    yLabel: r.radialSpec ? "Midpoint range r₀ (m)" : "Along-track offset x₀ (m)",
    yFormat: (x) => r.radialSpec ? x.toPrecision(3) : x.toFixed(2),
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
    yLabel: r.radialSpec ? "Radial acceleration a₀ (m/s²)" : "Perpendicular distance y₀ (m)",
    yFormat: (x) => r.radialSpec ? x.toPrecision(3) : x.toFixed(2),
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
  $("matchSummary").textContent =
    `${r.pulses} pulses · chirps ${r.start}–${r.start + r.pulses - 1} · verified peak ${db(r.best[3]).toFixed(2)} dB · ${r.radialSpec ? `r₀ ${x.toFixed(4)} m, v₀ ${y.toFixed(2)} m/s, a₀ ${v.toPrecision(5)} m/s² · FFT grid peak ${peak.toFixed(2)} dB` : `x₀ ${x.toFixed(4)} m, y₀ ${y.toFixed(4)} m, v ${v.toFixed(2)} m/s`} · ${noiseCount} noise-reference trains · ${(((r.pulses * meta.samples) / meta.parameters.fs) * 1e6).toFixed(2)} µs sampled / ${((r.pulses - 1) * r.period * 1e6 + (meta.samples / meta.parameters.fs) * 1e6).toFixed(2)} µs elapsed · ${r.seconds.toFixed(2)} s${r.radialSpec ? ` · grid ${nx} × ${nv} × ${ny} (r₀ × v₀ × a₀) · ${r.radialSpec[5]} correction FFTs · FFT ${r.radialSpec[3]} × ${r.radialSpec[4]}` : r.automatic ? ` · ${r.automatic[0]} automatic cells · FFT ${r.automatic[2]} × ${r.automatic[3]} · ${r.automatic[6]} shared FFT corrections · estimated cell phase loss ${r.loss}%` : " · direct grid"}${overlap ? " · selected train overlaps the mean-estimation interval" : ""}`;
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
      { name: "Re: measured residual", color: "#9badb9", x: xx, y: obs },
      {
        name: "Re: fitted template",
        color: "#f1cf74",
        x: xx,
        y: model,
        width: 1.7,
      },
    ],
    {
      x0: xx[0],
      x1: xx.at(-1),
      xLabel: "Time from train ramp start (µs; gaps are unsampled)",
      yLabel: "ADC counts",
    },
  );
  residualPlot.set(
    [
      {
        name: "Residual Re: data − fit",
        color: "#f1cf74",
        x: xx,
        y: Float64Array.from(obs, (value, i) => value - model[i]),
      },
      { name: "Residual Im", color: "#57d8c0", x: xx, y: residualImag },
    ],
    {
      x0: xx[0],
      x1: xx.at(-1),
      xLabel: "Time from train ramp start (µs; gaps are unsampled)",
      yLabel: "Residual (ADC counts)",
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
        color: "#f1cf74",
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
      yLabel: "Slant range (m)",
      decimals: 3,
    },
  );
  if (centeredZoom && num("zoom") > 0) {
    applyCenteredZoom(true);
  }
  if (r.scanResults) {
    const s = r.scanResults,
      starts = [],
      ys = [];
    for (let i = 0; i < s.length; i += 5) {
      starts.push(s[i]);
      ys.push(db(s[i + 1]));
    }
    scan.set(
      [
        {
          name: r.radialSpec ? "Time scan: verified radial FFT candidates" : "Time scan: maximum over x₀, y₀, v and I/Q orientation",
          color: "#57d8c0",
          x: starts,
          y: ys,
        },
      ],
      {
        x0: starts[0],
        x1: starts.at(-1) + 1,
        xLabel: "Train start chirp (click to select)",
        yLabel: "Energy / quiet noise (dB)",
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
          `<tr><td>${r.pulses}</td><td>${(((r.pulses * meta.samples) / meta.parameters.fs) * 1e6).toFixed(2)}</td><td>${db(r.best[3]).toFixed(2)}</td><td>${r.best[0].toFixed(4)}</td><td>${r.best[1].toFixed(4)}</td><td>${r.best[2].toFixed(2)}</td></tr>`,
      )
      .join("") +
    "</tbody></table>";
  const best = results.reduce((a, b) => (a.best[3] > b.best[3] ? a : b));
  showMatch(best);
  $("pulses").value = best.pulses;
  maps.forEach((m) => m.overlay({ bands: bands() }));
}
worker.onmessage = ({ data: m }) => {
  if (m.type === "loaded") {
    loaded = true;
    status(`RX${m.rx}: ${meta.total_chirps} chirps loaded; SHA-256 verified.`);
    setBusy(false);
  } else if (m.type === "background") {
    invalidate();
    power = m.power;
    const n = meta.total_chirps,
      mask = power.slice(2 * n),
      missing = [...mask].filter((x) => x === 0).length;
    $("padding").textContent =
      `${missing} / ${n} chirps flagged for zero runs ≥8 samples or nonfinite values. These are displayed in grey and excluded from noise estimates and searches. Original source reports ${meta.zero_padded_samples.toLocaleString()} padded samples across all receivers.`;
    $("bgNote").textContent =
      `Complex mean uses ${m.count} intact chirps from [${num("bgStart")}, ${num("bgStop")}). This same mean is removed from quiet and disturbed data.`;
    if (!busy) status("Complex background mean updated; plots refreshed.");
    applyTiming();
    view();
    updateSelection();
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
  } else if (m.type === "match") {
    showMatch(m.result);
  } else if (m.type === "complete") {
    setBusy(false);
    compareResults = m.results;
    if (m.results.length > 1) comparison(m.results);
    $("progress").value = 1;
    status(
      "Search complete. Peak values are referenced to the selected quiet background.",
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
$("background").onclick = safe(() => {
  assertRange(num("bgStart"), num("bgStop"));
  worker.postMessage({
    type: "background",
    start: num("bgStart"),
    stop: num("bgStop"),
  });
});
$("view").onclick = safe(customView);
$("zoom").oninput = () => {
  clearTimeout(zoomTimer);
  zoomTimer = setTimeout(safe(applyCenteredZoom), 90);
};
$("zoom").onchange = safe(() => {
  clearTimeout(zoomTimer);
  applyCenteredZoom();
});
$("timing").onclick = safe(applyTiming);
$("component").onchange = safe(view);
$("fftSub").onchange = safe(view);
for (const id of ["rawScale", "subScale", "span"])
  $(id).onchange = safe(drawView);
$("rx").onchange = () => {
  invalidate();
  setBusy(true);
  loaded = false;
  state.rx = num("rx");
  worker.postMessage({
    type: "rx",
    rx: state.rx,
    bgStart: num("bgStart"),
    bgStop: num("bgStop"),
  });
};
$("chirp").onchange = () => selectChirp(num("chirp"));
$("chirpSlider").oninput = () => selectChirp(num("chirpSlider"));
$("pulses").onchange = () => {
  invalidate();
  selectChirp(state.chirp);
};
$("traceSub").onchange = updateSelection;
$("full").onclick = safe(() => {
  clearTimeout(zoomTimer);
  $("zoom").value = 0;
  for(const plot of [...maps,vx,vy])plot.resetZoom();
  applyCenteredZoom();
});
$("frame").onclick = safe(() => {
  $("viewStart").value =
    Math.floor(state.chirp / meta.chirps_per_frame) * meta.chirps_per_frame;
  $("viewStop").value = num("viewStart") + meta.chirps_per_frame - 1;
  customView();
});
$("around").onclick = safe(() => {
  $("viewStart").value = Math.max(0, state.chirp - 125);
  $("viewStop").value = Math.min(meta.total_chirps - 1, state.chirp + 249);
  customView();
});
$("zoomBackground").onclick = safe(() => {
  $("viewStart").value = num("bgStart");
  $("viewStop").value = num("bgStop") - 1;
  customView();
});
$("search").onclick = () => run();
$("compare").onclick = () => run(true);
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
  "noiseStart",
  "noiseStop",
  "receiver",
  "phaseLoss",
])
  $(id).addEventListener("change", invalidate);
const savedBounds={fft:[.05,.7,0,1000000,-400,400],direct:[-.6,.6,.05,.3,200,500]};
let activeMethod='fft';
function methodUI() {
 const radial=$("algorithm").value==='fft';
 ['xN','yN','vN'].forEach((id,i)=>$(id).value=radial?'':[31,11,41][i]);
 $("axisR").textContent=radial?'r₀ (m)':'x₀ (m)';$("axisA").textContent=radial?'a₀ (m/s²)':'y₀ (m)';$("axisV").textContent=radial?'v₀ (m/s)':'v (m/s)';
 $("yMin").min=radial?'':'.001';$("yMax").min=radial?'':'.001';$("yMin").step=$("yMax").step=radial?'1000':'.01';
 $("mapRTitle").textContent=radial?'v₀ × r₀ · maximum over a₀':'v × x₀ · maximum over y₀';$("mapATitle").textContent=radial?'v₀ × a₀ · maximum over r₀':'v × y₀ · maximum over x₀';
 $("compare").textContent=radial?'Compare 2 / 4 / 8 / 16':'Compare 1 / 2 / 4 / 8 / 16';
 $("pulses").options[0].disabled=radial;if(radial && num('pulses')===1) $("pulses").value=8;
 $("modelNote").textContent=radial?'Polynomial radial range: r(t) = r₀ + v₀ h + ½a₀h². t₀ is the sampled train midpoint; h = t − t₀. Positive v₀ is receding. a₀ is radial acceleration (m/s²); the carrier Doppler chirp rate is approximately −2f₀a₀/c Hz/s. These are independent search coordinates.':'Radar at (x₀, y₀), projectile at (v t, 0). x₀ is along-track offset at this train’s start; y₀ is perpendicular distance. R(t) = √((v t − x₀)² + y₀²).';
 $("gridNote").textContent=radial?'Enter r₀, v₀ and a₀ bounds. Acceleration spacing and shared velocity corrections follow a polynomial phase bound over the full coherent train. Range/velocity spacing comes from padded FFT bins. The grid scores approximate the phase-only quadratic matched filter; strongest candidates are verified directly. This model omits cubic and higher range terms. Optional HPFs apply to candidate verification; the FFT grid is phase-only.':'Enter x₀, y₀ and signed along-track velocity bounds and grid counts for the full retarded Eq. 25 reference search.';
 $("mapNote").textContent=radial?'MAX projections of the complete r₀ / v₀ / a₀ FFT grid, including both I/Q orientations. Gray cells have nonpositive range over the train. Crosses mark the FFT grid peak. Strong candidates are checked with the exact quadratic template for the reported fit. FFT-bin rounding and correction sharing make grid scores approximate. A bank peak alone does not identify a projectile.':'MAX projections over the full configured geometry grid and both I/Q orientations. Crosses mark the global best fit. Mirroring (x₀,v) gives the same monostatic range history. A bank peak alone does not identify a projectile.';
}
$("algorithm").onchange = () => {
 savedBounds[activeMethod]=['xMin','xMax','yMin','yMax','vMin','vMax'].map(num);activeMethod=$("algorithm").value;
 ['xMin','xMax','yMin','yMax','vMin','vMax'].forEach((id,i)=>$(id).value=savedBounds[activeMethod][i]);
 invalidate();methodUI();setBusy(false);
};
methodUI();
try {
  setBusy(false);
  const response = await fetch("datasets/test63/metadata.json", {
    cache: "no-store",
  });
  if (!response.ok) throw Error("Recording metadata unavailable");
  meta = await response.json();
  clockBase = Date.parse(meta.recorded_clock + "Z");
  $("record").textContent =
    `Test 63 · ${meta.diameter_mm} mm ${meta.material} ball · ${meta.parameters.speed.toFixed(2)} m/s · ${(meta.parameters.f_start / 1e9).toFixed(1)} GHz · ${meta.frames} frames × ${meta.chirps_per_frame} chirps × ${meta.samples} samples × ${meta.receivers} RX`;
  $("viewStop").value = meta.total_chirps - 1;
  $("chirpSlider").max = meta.total_chirps - 1;
  $("framePeriod").value = (
    (meta.parameters.T_adc +
      meta.samples / meta.parameters.fs +
      meta.parameters.T_idle +
      meta.ramp_tail_us_assumption * 1e-6) *
    meta.chirps_per_frame *
    1e3
  ).toFixed(5);
  $("fftNote").textContent =
    `Complex Hann-window FFT: ${meta.samples} acquired samples, ${2 ** Math.ceil(Math.log2(meta.samples))} FFT points. Fourier resolution ${(meta.parameters.fs / meta.samples / 1000).toFixed(2)} kHz; zero-padded bin spacing ${(meta.parameters.fs / 2 ** Math.ceil(Math.log2(meta.samples)) / 1000).toFixed(2)} kHz. Both spectrum views share the original recording's peak reference; frequency includes range beat and Doppler.`;
  worker.postMessage({
    type: "init",
    meta,
    bgStart: num("bgStart"),
    bgStop: num("bgStop"),
  });
} catch (e) {
  status(e.message, true);
  console.error(e);
}
