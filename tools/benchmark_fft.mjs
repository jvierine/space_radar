// Compare the shipped Wasm search paths on identical data and physical bounds.
import fs from "node:fs";
const input = JSON.parse(fs.readFileSync(0,"utf8"));
const {instance} = await WebAssembly.instantiate(fs.readFileSync("web/lab/core.wasm"),{});
const w=instance.exports;
const read=()=>Array.from(new Float32Array(w.memory.buffer,w.result_ptr(),w.result_len()));
const p=input.p, data=new Float32Array(input.data),ptr=w.allocate(data.length);
new Float32Array(w.memory.buffer,ptr,data.length).set(data);
w.load(ptr,p.rows,p.n,p.frame,p.fs,p.adc,p.period,p.f0,p.slope);w.release(ptr,data.length);
w.background(0,125);
const reports=[];
for(const pulses of input.pulses ?? [4,8,16]) {
  w.prepare(750,pulses,125,625,input.receiver ? 1 : 0);
  const bounds=input.bounds ?? [-.6,.6,.05,.3,200,500];
  let time=performance.now();
  w.auto_begin(...bounds,input.loss ?? .05,150000);
  let cells=0,steps=0;
  while(!cells && steps++<10000) {
    cells=w.auto_build_batch(128);
    if(cells<0) throw Error(`Auto grid failed: ${cells}`);
  }
  if(!cells) throw Error("Automatic grid never completed");
  w.auto_info(); const specification=read();
  const build_seconds=(performance.now()-time)/1000;
  console.error(`${pulses} pulses: ${cells} cells, built in ${build_seconds.toFixed(3)} s`);
  w.auto_cells(); const boxes=read();
  w.auto_nodes(); const initial_nodes=read();
  time=performance.now();
  for(let first=0;first<cells;first+=1) w.search_batch(first,1);
  const fft_seconds=(performance.now()-time)/1000;
  time=performance.now(); w.auto_refine(32);
  const refine_seconds=(performance.now()-time)/1000;
  w.matches();const packed=read(),best=packed.slice(cells),scores=packed.slice(0,cells);
  w.auto_nodes();const nodes=read();
  w.fitted();const fitted=read();
  w.train_values();const observed=read();
  time=performance.now();
  const count=w.grid(bounds[0],bounds[1],31,bounds[2],bounds[3],11,bounds[4],bounds[5],41);
  w.search_batch(0,count); w.matches();const direct=read().slice(count);
  const direct_seconds=(performance.now()-time)/1000;
  console.error(`${pulses} pulses: FFT ${fft_seconds.toFixed(3)} s + refinement ${refine_seconds.toFixed(3)} s; direct ${direct_seconds.toFixed(3)} s`);
  // Exact physical evaluation at the returned point, in both orientations.
  w.grid(best[0],best[0],1,best[1],best[1],1,best[2],best[2],1);w.search_batch(0,1);w.matches();
  const verified=read().slice(1);
  reports.push({pulses,specification,build_seconds,fft_seconds,refine_seconds,direct_seconds,
    best,direct,verified,boxes,nodes,initial_nodes,scores,fitted,observed});
}
console.log(JSON.stringify(reports));
