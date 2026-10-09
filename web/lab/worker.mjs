import {projectReceivers,phaseSearch} from './beamforming.mjs?v=20261009rxnoise10';
import {commonValid,prepareReceiverTrains} from './receiver-trains.mjs?v=20261009rxnoise10';
import {runRadial} from './radial-backend.mjs?v=20261009rxnoise10';
let wasm,
  meta,
  data,
  rx = 0,
  generation = 0,
  bgRange,
  noiseRange;
let searchFinished=Promise.resolve();
const streams=new Map();
const send = (type, rest = {}) => postMessage({ type, ...rest });
const copy = () =>
  new Float32Array(
    wasm.memory.buffer,
    wasm.result_ptr(),
    wasm.result_len(),
  ).slice();
const yieldUI = () => new Promise((r) => setTimeout(r, 0));
async function channel(receiver) {
  if(streams.has(receiver))return streams.get(receiver);
  const spec=meta.transport[receiver];
  const response=await fetch('datasets/test63/'+spec.file);
  if(!response.ok)throw Error('Cannot load receiver data');
  const bytes=await response.arrayBuffer();
  if(bytes.byteLength!==spec.bytes)throw Error('Receiver file length mismatch');
  const sha=[...new Uint8Array(await crypto.subtle.digest('SHA-256',bytes))].map(x=>x.toString(16).padStart(2,'0')).join('');
  if(sha!==spec.sha256)throw Error('Receiver SHA-256 mismatch');
  const result=new Float32Array(bytes);streams.set(receiver,result);return result;
}
async function load(receiver) {
  rx=receiver;data=await channel(rx);
  const ptr = wasm.allocate(data.length);
  new Float32Array(wasm.memory.buffer, ptr, data.length).set(data);
  const p = meta.parameters;
  wasm.load(
    ptr,
    meta.total_chirps,
    meta.samples,
    meta.chirps_per_frame,
    p.fs,
    p.T_adc,
    p.T_adc +
      meta.samples / p.fs +
      p.T_idle +
      meta.ramp_tail_us_assumption * 1e-6,
    p.f_start,
    p.freq_slope,
  );
  wasm.release(ptr, data.length);
  send("loaded", { rx, verified_sha256: meta.transport[rx].sha256 });
}
function setBackground(a, b) {
  bgRange = [a, b];
  const count = wasm.background(a, b);
  if (!count) throw Error("Quiet mean interval contains no valid chirps");
  wasm.powers();
  const power = copy();
  send("background", {
    count,
    power,
    fft_n: 2 ** Math.ceil(Math.log2(meta.samples)),
  });
}
async function view(start, stop, component, fftSub) {
  const images = [];
  for (const mode of [component, component + 1, fftSub ? 5 : 4]) {
    const height = wasm.image(start, stop, mode);
    images.push({ mode, height, values: copy() });
  }
  wasm.image(0,meta.total_chirps,4);const ref=copy();
  let peak=0;for(const v of ref)if(Number.isFinite(v)&&v>peak)peak=v;
  send("view", { start, stop, images, fft_peak: peak });
}
function trace(chirp) {
  const valid = wasm.trace(chirp);
  send("trace", { chirp, valid, values: copy() });
}
function setup(job, pulses, start) {
  noiseRange = [job.noiseStart, job.noiseStop];
  const n = wasm.prepare_power(start, pulses, ...noiseRange, +job.receiver);
  if (n === -1)
    throw Error(
      "This train crosses a frame boundary or contains padded samples. Select an intact train.",
    );
  if (n === -2)
    throw Error(
      "Blue background needs one intact chirp for full-bandwidth noise power.",
    );
  const g = job.grid;
  if(job.algorithm === "fft") return 0;
  if (![g.xN,g.yN,g.vN].every(n=>Number.isInteger(n)&&n>0) || g.xN*g.yN*g.vN>150000) throw Error("Direct reference grid must contain 1–150,000 templates.");
  return wasm.grid(
    g.xMin,
    g.xMax,
    g.xN,
    g.yMin,
    g.yMax,
    g.yN,
    g.vMin,
    g.vMax,
    g.vN,
  );
}
async function bank(job, pulses, start, token, phase = "search", scan = false) {
  const started = performance.now();
  let total = setup(job, pulses, start);
  const receivers=await Promise.all(meta.transport.map((_,i)=>channel(i)));
  if(token!==generation)return null;
  const trains=prepareReceiverTrains(receivers,{samples:meta.samples,rows:meta.total_chirps,perFrame:meta.chirps_per_frame,start,pulses,bgStart:bgRange[0],bgStop:bgRange[1],noiseStart:job.noiseStart,noiseStop:job.noiseStop,fullBandwidth:true});
  const ptr=wasm.allocate(trains.data.length);
  try {
    new Float32Array(wasm.memory.buffer,ptr,trains.data.length).set(trains.data);
    if(wasm.search_receivers(ptr,trains.data.length,4,trains.quietCount,rx)<0)throw Error('Invalid four-receiver trains');
    const np=wasm.allocate(4);new Float32Array(wasm.memory.buffer,np,4).set(trains.noisePower);
    try{if(wasm.search_noise_power(np,4)<0)throw Error('Invalid raw background noise power');}finally{wasm.release(np,4);}
  } finally {wasm.release(ptr,trains.data.length);}

  let automatic = null, cells = null;
  let radialSpec=null,compute={backend:"CPU (Rust/Wasm)",fallbackReason:null};
  if (job.algorithm === "fft") {
    const g=job.grid;
    total=wasm.radial_begin(g.xMin,g.xMax,g.vMin,g.vMax,g.yMin,g.yMax,job.loss/100,32000000);
    if(total<0) throw Error("Invalid radial bounds or more than 32 million r₀ / v₀ / a₀ grid points. Narrow the bounds or increase the phase tolerance.");
    wasm.radial_info();radialSpec=[...copy()];
  }
  let windows = 0;
  if (scan && !radialSpec) {
    if ((job.scanStop - job.scanStart) / job.scanStride > 512)
      throw Error(
        "Limit each time scan to 512 starts; increase stride or shorten the interval.",
      );
    windows = wasm.scan_range(job.scanStart, job.scanStop, job.scanStride);
    if (!windows)
      throw Error("No intact independent trains in the scan interval.");
  }
  if(radialSpec){
    const g=job.grid;
    try{compute=await runRadial({w:wasm,info:radialSpec,backend:job.computeBackend??'auto',
      begin:()=>wasm.radial_begin(g.xMin,g.xMax,g.vMin,g.vMax,g.yMin,g.yMax,job.loss/100,32000000),
      cancelled:()=>token!==generation,yieldUI,
      progress:(fraction,backend)=>send('progress',{fraction,pulses,phase:`${phase} · ${backend} · ${total} corrections`,windows})});
    }catch(error){if(token!==generation)return null;throw error;}
  }
  const batch = 32;
  for (let first = 0; first < (radialSpec ? 0 : total); first += batch) {
    if (token !== generation) return null;
    if(radialSpec) wasm.radial_batch(first,batch); else wasm.search_batch(first, batch);
    send("progress", {
      fraction: Math.min((first + batch) / total, 1),
      pulses,
      phase: `${phase} · ${total} ${radialSpec ? "radial acceleration / velocity corrections" : "templates"}`,
      windows,
    });
    await yieldUI();
  }
  if (radialSpec) {
    send("progress",{fraction:1,pulses,phase:"Verifying radial FFT candidates with the quadratic template"});
    await yieldUI();if(token!==generation)return null;if(wasm.radial_refine(128)<0) throw Error("No valid positive-range quadratic templates in these bounds.");
  }
  if (automatic) {
    send("progress", {fraction:1,pulses,phase:"Refining FFT peaks with the physical template",windows});
    await yieldUI();
    if (token !== generation) return null;
    wasm.auto_refine(32);
  }
  let nodes=null;
  if(automatic){wasm.auto_nodes();nodes=copy();}
  const cubeLength = wasm.matches(),
    packed = copy();
  wasm.fitted();
  const fit = copy();
  wasm.train_values();
  const observed = copy();
  let scanResults = null;
  if (scan) {
    wasm.scan_results();
    scanResults = copy();
  }
  return {
    pulses,
    start,
    cube: packed.slice(0, cubeLength),
    best: [...packed.slice(cubeLength)],
    fit,
    observed,
    scanResults,
    grid: radialSpec ? {...job.grid,xN:radialSpec[0],yN:radialSpec[1],vN:radialSpec[2]} : job.grid,
    radialSpec,
    model: radialSpec ? "radial-quadratic" : "geometry",
    midpoint: meta.parameters.T_adc+(meta.samples-1)/(2*meta.parameters.fs)+(pulses-1)*job.period/2,
    noiseRange,
    bgRange,
    receiver: job.receiver,
    searchReceivers: 4,
    referenceStarts: trains.referenceStarts,
    noisePower:trains.noisePower,noiseSamples:trains.noiseSamples,noiseModel:"raw full-bandwidth sample power",
    rx,
    period: job.period,
    algorithm: job.algorithm ?? "direct",
    ...compute,
    automatic,
    cells,
    nodes,
    seconds: (performance.now()-started)/1000,
    loss: job.loss,
  };
}
onmessage = async ({ data: msg }) => {
  try {
    if (msg.type === "init") {
      meta = msg.meta;
      const result = await WebAssembly.instantiateStreaming(
        fetch("core.wasm?v=20261009rxnoise10"),
        {},
      );
      wasm = result.instance.exports;
      await load(msg.rx ?? 0);
      setBackground(msg.bgStart, msg.bgStop);
    } else if (msg.type === "rx") {
      generation++;
      await load(msg.rx);
      setBackground(msg.bgStart, msg.bgStop);
    } else if (msg.type === "background") {
      generation++;
      setBackground(msg.start, msg.stop);
    } else if (msg.type === "period") {
      wasm.set_period(msg.period);
    } else if (msg.type === "view") {
      await view(msg.start, msg.stop, msg.component, msg.fftSub);
    } else if (msg.type === "trace") {
      trace(msg.chirp);
    } else if (msg.type === "cancel") {
      generation++;
      send("cancelled");
    } else if (msg.type === "search") {
      const token = ++generation;
      const previous=searchFinished;
      let release;
      searchFinished=new Promise(resolve=>{release=resolve;});
      await previous;
      try {
      if(token!==generation)return;
      const results = [];
      for (const n of msg.compare ? (msg.algorithm==="fft" ? [2,4,8,16] : [1,2,4,8,16]) : [msg.pulses]) {
        if(msg.scan) {
          let winner=null,scan=[],overlapSkipped=0,invalidSkipped=0;
          const receivers=await Promise.all(meta.transport.map((_,i)=>channel(i)));
          const valid=commonValid(receivers,meta.samples,meta.total_chirps);
          send('scan-start',{start:msg.scanStart,stop:msg.scanStop});
          for(let start=msg.scanStart;start+n<=msg.scanStop;start++) {
            if(token!==generation)return;
            if(start<msg.noiseStop && start+n>msg.noiseStart || start<bgRange[1] && start+n>bgRange[0]){overlapSkipped++;continue;}
            if(Math.floor(start/meta.chirps_per_frame)!==Math.floor((start+n-1)/meta.chirps_per_frame)||!valid.slice(start,start+n).every(Boolean)){invalidSkipped++;continue;}
            const prepared=wasm.prepare_power(start,n,msg.noiseStart,msg.noiseStop,+msg.receiver);
            if(prepared===-2)throw Error('Blue background needs one intact chirp for full-bandwidth noise power.');
            if(prepared<0){invalidSkipped++;continue;}
            const candidate=await bank({...msg,scan:false},n,start,token,'time scan');
            if(!candidate)return;
            const p=meta.parameters;
            (candidate.model==='radial-quadratic'?wasm.radial_template:wasm.template_values)(...candidate.best.slice(0,3),n,+msg.receiver);
            const projections=projectReceivers(receivers,copy(),{samples:meta.samples,rows:meta.total_chirps,perFrame:meta.chirps_per_frame,start,pulses:n,bgStart:bgRange[0],bgStop:bgRange[1],noiseStart:msg.noiseStart,noiseStop:msg.noiseStop,conjugated:candidate.best[4]>0,fullBandwidth:true});
            candidate.beam=phaseSearch(projections.event,projections.noise,msg.beamSteps??10,projections.noiseCovariance);
            // RCS sampled coherent time is unchanged by receiver phases.
            (candidate.model==='radial-quadratic'?wasm.radial_template:wasm.template_values)(...candidate.best.slice(0,3),n,+msg.receiver);
            const q=copy();let sum=0,energy=0;
            for(let i=0;i<q.length;i+=2){const a=Math.hypot(q[i],q[i+1]);sum+=a;energy+=a*a;}
            candidate.beam.effectiveTime=sum*sum/(p.fs*energy);
            candidate.beam.rawNoisePower=projections.rawNoisePower;candidate.beam.noiseSamples=projections.noiseSamples;
            candidate.beam.referenceStarts=projections.referenceStarts;candidate.beam.meanCount=projections.meanCount;
            candidate.scanPoint=true;send('match',{result:candidate});
            scan.push(start,candidate.best[3],...candidate.best.slice(0,3));
            send('scan-progress',{fraction:(start+n-msg.scanStart)/(msg.scanStop-msg.scanStart),start});
            await yieldUI();if(token!==generation)return;
            if(!winner || candidate.best[3]>winner.best[3])winner=candidate;
          }
          if(!winner)throw Error(`No usable ${n}-chirp trains in the purple analysis window: ${overlapSkipped} overlap the background; ${invalidSkipped} contain padding or cross a frame. Move or enlarge the purple window.`);
          winner.scanResults=Float32Array.from(scan);winner.scanPoint=false;results.push(winner);send('match',{result:winner});continue;
        }
        let result = await bank(
          msg,
          n,
          msg.start,
          token,
          "selected train",
          msg.scan,
        );
        if (!result) return;
        if (msg.scan && result.scanResults?.length) {
          const scan = result.scanResults;
          let best = 0;
          for (let i = 5; i < scan.length; i += 5)
            if (scan[i + 1] > scan[best + 1]) best = i;
          const peakStart = scan[best];
          result = await bank(
            msg,
            n,
            peakStart,
            token,
            "map at time-scan peak",
          );
          if (!result) return;
          result.scanResults = scan;
        }
        results.push(result);
        send("match", { result, interim: !!msg.compare });
      }
      send("complete", { results });
      } finally {release();}
    }
  } catch (error) {
    send("error", { message: error.message, stack: error.stack });
  }
};
