let wasm,
  meta,
  data,
  rx = 0,
  generation = 0,
  bgRange,
  noiseRange;
const send = (type, rest = {}) => postMessage({ type, ...rest });
const copy = () =>
  new Float32Array(
    wasm.memory.buffer,
    wasm.result_ptr(),
    wasm.result_len(),
  ).slice();
const yieldUI = () => new Promise((r) => setTimeout(r, 0));
async function load(receiver) {
  rx = receiver;
  const spec = meta.transport[rx];
  const response = await fetch("datasets/test63/" + spec.file);
  if (!response.ok) throw Error("Cannot load receiver data");
  const bytes = await response.arrayBuffer();
  if (bytes.byteLength !== spec.bytes)
    throw Error("Receiver file length mismatch");
  const sha = [...new Uint8Array(await crypto.subtle.digest("SHA-256", bytes))]
    .map((x) => x.toString(16).padStart(2, "0"))
    .join("");
  if (sha !== spec.sha256) throw Error("Receiver SHA-256 mismatch");
  data = new Float32Array(bytes);
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
  send("loaded", { rx, verified_sha256: sha });
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
function view(start, stop, component, fftSub) {
  const images = [];
  for (const mode of [component, component + 1, fftSub ? 5 : 4]) {
    const height = wasm.image(start, stop, mode);
    images.push({ mode, height, values: copy() });
  }
  wasm.image(0, meta.total_chirps, 4);
  const ref = copy();
  let peak = 0;
  for (const v of ref) if (Number.isFinite(v) && v > peak) peak = v;
  send("view", { start, stop, images, fft_peak: peak });
}
function trace(chirp) {
  const valid = wasm.trace(chirp);
  send("trace", { chirp, valid, values: copy() });
}
function setup(job, pulses, start) {
  noiseRange = [job.noiseStart, job.noiseStop];
  const n = wasm.prepare(start, pulses, ...noiseRange, +job.receiver);
  if (n === -1)
    throw Error(
      "This train crosses a frame boundary or contains padded samples. Select an intact train.",
    );
  if (n === -2)
    throw Error(
      "Need at least eight intact, independent noise-reference trains for this pulse count. Enlarge the noise interval.",
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
  if (pulses === 1 && job.algorithm === "fft") job = {...job, algorithm:"direct"};
  const started = performance.now();
  let total = setup(job, pulses, start);
  let automatic = null, cells = null;
  if (job.algorithm === "fft") {
    const g = job.grid;
    if (wasm.auto_begin(g.xMin,g.xMax,g.yMin,g.yMax,g.vMin,g.vMax,job.loss/100,150000) < 0)
      throw Error("Invalid bounds or correction-phase loss limit.");
    for (;;) {
      if (token !== generation) return null;
      const done = wasm.auto_build_batch(128);
      if (done === -1) throw Error("The requested correction grid exceeds 150,000 cells. Narrow the bounds, shorten the coherent train, or increase the phase-loss allowance.");
      if (done < 0) throw Error("Cannot bound corrections over this domain. Use positive y₀ and a finite velocity interval.");
      wasm.auto_info(); automatic = [...copy()];
      send("progress", {fraction:0, pulses,phase:`Building automatic grid · ${automatic[0]} accepted / ${automatic[1]} pending`});
      if (done > 0) {total=done;break;}
      await yieldUI();
    }
    wasm.auto_cells(); cells=copy();
  }
  let windows = 0;
  if (scan) {
    if ((job.scanStop - job.scanStart) / job.scanStride > 512)
      throw Error(
        "Limit each time scan to 512 starts; increase stride or shorten the interval.",
      );
    windows = wasm.scan_range(job.scanStart, job.scanStop, job.scanStride);
    if (!windows)
      throw Error("No intact independent trains in the scan interval.");
  }
  const batch = job.algorithm === "fft" ? 8 : 32;
  for (let first = 0; first < total; first += batch) {
    if (token !== generation) return null;
    wasm.search_batch(first, batch);
    send("progress", {
      fraction: Math.min((first + batch) / total, 1),
      pulses,
      phase: `${phase} · ${total} ${automatic ? "midpoint correction cells" : "templates"}`,
      windows,
    });
    await yieldUI();
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
    grid: job.grid,
    noiseRange,
    bgRange,
    receiver: job.receiver,
    rx,
    period: job.period,
    algorithm: job.algorithm ?? "direct",
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
        fetch("core.wasm?v=20261009midpoint1"),
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
      view(msg.start, msg.stop, msg.component, msg.fftSub);
    } else if (msg.type === "trace") {
      trace(msg.chirp);
    } else if (msg.type === "cancel") {
      generation++;
      send("cancelled");
    } else if (msg.type === "search") {
      const token = ++generation;
      const results = [];
      for (const n of msg.compare ? [1, 2, 4, 8, 16] : [msg.pulses]) {
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
    }
  } catch (error) {
    send("error", { message: error.message, stack: error.stack });
  }
};
