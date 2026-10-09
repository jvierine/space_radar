import { Heatmap, LinePlot, db } from "./plots.mjs?v=20261009b";
const $ = (id) => document.getElementById(id),
  num = (id) => Number($(id).value),
  worker = new Worker("worker.mjs?v=20261009b", { type: "module" });
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
  clockBase;
const state = { start: 0, stop: 6250, chirp: 625, pulses: 4, rx: 0 };
function status(text, error = false) {
  $("status").textContent = text;
  $("status").classList.toggle("error", error);
}
function setBusy(value) {
  busy = value;
  for (const el of document.querySelectorAll(
    ".settings input,.settings select,.settings button,#matchSettings input,#matchSettings select,#search,#compare,#scan",
  ))
    el.disabled = value || !loaded;
  $("cancel").disabled = !value;
}
function assertRange(a, b) {
  if (
    !Number.isInteger(a) ||
    !Number.isInteger(b) ||
    a < 0 ||
    b > meta.total_chirps ||
    a >= b
  )
    throw Error("Enter a valid, nonempty chirp interval.");
}
function clock(chirp) {
  const k = Math.round(chirp),
    frame = Math.floor(k / meta.chirps_per_frame),
    within = k % meta.chirps_per_frame;
  const elapsed = frame * framePeriod + within * period + meta.parameters.T_adc;
  const ms = clockBase + elapsed * 1000,
    d = new Date(Math.floor(ms));
  return `${d.toISOString().slice(11, 19)}.${Math.floor((ms % 1000) * 1000)
    .toString()
    .padStart(6, "0")}`;
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
  scan = new LinePlot("scanPlot", (k) => selectChirp(Math.round(k)));
function common() {
  return {
    x0: state.start,
    x1: state.stop,
    y0: meta.parameters.T_adc * 1e6,
    y1: (meta.parameters.T_adc + (meta.samples - 1) / meta.parameters.fs) * 1e6,
    xFormat: (x) => Math.min(state.stop - 1, Math.round(x)).toString(),
    topFormat: (x) => clock(Math.min(state.stop - 1, Math.round(x))),
    topLabel: "Wall clock (reconstructed; timezone unspecified)",
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
    b = num("viewStop");
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
function invalidate() {
  result = undefined;
  vx.root.style.display = "none";
  vy.root.style.display = "none";
  fit.root.style.display = "none";
  scan.root.style.display = "none";
  $("compareResults").innerHTML = "";
  $("matchSummary").textContent =
    "Settings changed: run a search to update the filter maps.";
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
  const desc = `Global ${msg.chirp} · frame ${Math.floor(msg.chirp / meta.chirps_per_frame)}, chirp ${msg.chirp % meta.chirps_per_frame} · reconstructed clock ${clock(msg.chirp)}${msg.valid ? "" : " · FLAGGED PADDING: excluded from filters/spectra"}`;
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
      !Number.isInteger(g[a + "N"]) ||
      g[a + "N"] < 1
    )
      throw Error("Invalid trajectory grid.");
  if (g.yMin <= 0 || Math.max(Math.abs(g.vMin), Math.abs(g.vMax)) > 100000)
    throw Error("Use positive perpendicular distance and |v| ≤ 100 km/s.");
  const cells = g.xN * g.yN * g.vN;
  if (cells > 150000)
    throw Error(
      "Limit the bank to 150,000 templates; use a coarse search then refine.",
    );
  return g;
}
function run(compare = false, timeScan = false) {
  if (busy) return;
  try {
    applyTiming();
    assertRange(num("bgStart"), num("bgStop"));
    worker.postMessage({
      type: "background",
      start: num("bgStart"),
      stop: num("bgStop"),
    });
    const g = grid();
    assertRange(num("noiseStart"), num("noiseStop"));
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
      compare,
      scan: timeScan,
      scanStart: num("scanStart"),
      scanStop: num("scanStop"),
      scanStride: num("scanStride"),
    };
    if (timeScan) {
      assertRange(job.scanStart, job.scanStop);
      if (!Number.isInteger(job.scanStride) || job.scanStride < 1)
        throw Error("Start stride must be a positive integer.");
      if ((job.scanStop - job.scanStart) / job.scanStride > 512)
        throw Error("Use at most 512 time starts per scan; increase stride.");
    }
    setBusy(true);
    compareResults = [];
    $("compareResults").innerHTML = "";
    $("progress").value = 0;
    status(`Searching ${g.xN * g.yN * g.vN} trajectory templates…`);
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
  result = r;
  state.pulses = r.pulses;
  vx.root.style.display = "block";
  vy.root.style.display = "block";
  maps.forEach((m) => m.overlay({ bands: bands() }));
  const g = r.grid,
    nv = g.vN,
    nx = g.xN,
    ny = g.yN,
    a = new Float32Array(nx * nv).fill(-Infinity),
    b = new Float32Array(ny * nv).fill(-Infinity);
  for (let ix = 0; ix < nx; ix++)
    for (let iy = 0; iy < ny; iy++)
      for (let iv = 0; iv < nv; iv++) {
        const score = db(r.cube[(ix * ny + iy) * nv + iv]);
        if (score > a[ix * nv + iv]) a[ix * nv + iv] = score;
        if (score > b[iy * nv + iv]) b[iy * nv + iv] = score;
      }
  const peak = db(r.best[3]),
    low = Math.min(0, peak - 35),
    high = Math.max(10, peak);
  const opts = {
    x0: g.vMin,
    x1: g.vMax === g.vMin ? g.vMax + 1 : g.vMax,
    lo: low,
    hi: high,
    kind: 1,
    xLabel: "Along-track velocity v (m/s)",
    xFormat: (x) => x.toFixed(0),
    colorLabel: "Matched energy / quiet noise (dB)",
  };
  vx.set(a, nv, nx, {
    ...opts,
    y0: g.xMin,
    y1: g.xMax === g.xMin ? g.xMax + 0.01 : g.xMax,
    yLabel: "Along-track offset x₀ (m)",
    yFormat: (x) => x.toFixed(2),
  });
  vy.set(b, nv, ny, {
    ...opts,
    y0: g.yMin,
    y1: g.yMax === g.yMin ? g.yMax + 0.01 : g.yMax,
    yLabel: "Perpendicular distance y₀ (m)",
    yFormat: (x) => x.toFixed(2),
  });
  const [x, y, v, , , noiseCount] = r.best;
  const overlap =
    r.start < num("bgStop") && r.start + r.pulses > num("bgStart");
  $("matchSummary").textContent =
    `${r.pulses} pulses · chirps ${r.start}–${r.start + r.pulses - 1} · peak ${peak.toFixed(2)} dB · x₀ ${x.toFixed(4)} m, y₀ ${y.toFixed(4)} m, v ${v.toFixed(2)} m/s · ${noiseCount} noise-reference trains · ${(((r.pulses * meta.samples) / meta.parameters.fs) * 1e6).toFixed(2)} µs sampled / ${((r.pulses - 1) * r.period * 1e6 + (meta.samples / meta.parameters.fs) * 1e6).toFixed(2)} µs elapsed${overlap ? " · selected train overlaps the mean-estimation interval" : ""}`;
  const xx = [],
    obs = [],
    model = [];
  for (let k = 0; k < r.pulses; k++) {
    for (let j = 0; j < meta.samples; j++) {
      xx.push(
        (k * r.period + meta.parameters.T_adc + j / meta.parameters.fs) * 1e6,
      );
      obs.push(r.observed[2 * (k * meta.samples + j)]);
      model.push(r.fit[2 * (k * meta.samples + j)]);
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
          name: "Time scan: maximum over x₀, y₀, v and I/Q orientation",
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
  $("compareResults").innerHTML =
    "<table><thead><tr><th>Pulses</th><th>Sampled µs</th><th>Peak dB</th><th>x₀ (m)</th><th>y₀ (m)</th><th>v (m/s)</th></tr></thead><tbody>" +
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
$("view").onclick = safe(view);
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
  $("viewStart").value = 0;
  $("viewStop").value = meta.total_chirps;
  view();
});
$("frame").onclick = safe(() => {
  $("viewStart").value =
    Math.floor(state.chirp / meta.chirps_per_frame) * meta.chirps_per_frame;
  $("viewStop").value = num("viewStart") + meta.chirps_per_frame;
  view();
});
$("around").onclick = safe(() => {
  $("viewStart").value = Math.max(0, state.chirp - 125);
  $("viewStop").value = Math.min(meta.total_chirps, state.chirp + 250);
  view();
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
])
  $(id).addEventListener("change", invalidate);
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
  $("viewStop").value = meta.total_chirps;
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
