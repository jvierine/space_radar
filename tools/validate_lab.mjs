// Numerical harness for the exact compiled browser Wasm (not a JS substitute).
import fs from "node:fs";
const input = JSON.parse(fs.readFileSync(0, "utf8"));
const { instance } = await WebAssembly.instantiate(
    fs.readFileSync("web/lab/core.wasm"),
    {},
  ),
  w = instance.exports;
const read = () =>
  Array.from(new Float32Array(w.memory.buffer, w.result_ptr(), w.result_len()));
function load(values, p) {
  const a = new Float32Array(values),
    ptr = w.allocate(a.length);
  new Float32Array(w.memory.buffer, ptr, a.length).set(a);
  w.load(ptr, p.rows, p.n, p.frame, p.fs, p.adc, p.period, p.f0, p.slope);
  w.release(ptr, a.length);
}
load(input.data, input.p);
const baselineCount = w.background(0, 125);
w.powers();
const powers = read();
w.trace(750);
const trace = read();
const templates = [];
for (const c of input.cases) {
  w.template_values(c.x, c.y, c.v, c.pulses, c.receiver);
  templates.push(read());
}
const scores = [];
for (const n of [1, 2, 4, 8, 16]) {
  const refs = w.prepare(750, n, 125, 625, 0);
  if (refs < 8) throw Error("Insufficient noise reference");
  const cells = w.grid(0.29, 0.31, 3, 0.09, 0.11, 3, 331.51, 351.51, 3);
  w.search_batch(0, cells);
  w.matches();
  scores.push({ n, values: read() });
}
const paddingRejected = w.prepare(880, 4, 125, 625, 0),
  crossFrameRejected = w.prepare(124, 4, 125, 625, 0);
w.prepare(750, 4, 125, 625, 0);
w.grid(0.3, 0.3, 1, 0.1, 0.1, 1, 341.51, 341.51, 1);
const scanCount = w.scan_range(0, 1000, 1);
w.search_batch(0, 1);
w.scan_results();
const scanned = read();
console.log(
  JSON.stringify({
    baselineCount,
    powers,
    trace,
    templates,
    scores,
    paddingRejected,
    crossFrameRejected,
    scanCount,
    scanned,
  }),
);
