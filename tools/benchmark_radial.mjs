import fs from 'node:fs';
const job=JSON.parse(fs.readFileSync(0,'utf8'));
const w=(await WebAssembly.instantiate(fs.readFileSync('web/lab/core.wasm'),{})).instance.exports;
const read=()=>Array.from(new Float32Array(w.memory.buffer,w.result_ptr(),w.result_len()));
const p=job.p,data=Float32Array.from(job.data),ptr=w.allocate(data.length);
new Float32Array(w.memory.buffer,ptr,data.length).set(data);
w.load(ptr,p.rows,p.n,p.frame,p.fs,p.adc,p.period,p.f0,p.slope);w.release(ptr,data.length);w.background(0,125);
const prepared=w.prepare(750,job.pulses,125,625,job.receiver??0);if(prepared<0)throw Error('prepare '+prepared);
const t=performance.now();const count=w.radial_begin(...job.bounds,job.loss??.05,job.cap??32000000);if(count<0)throw Error('radial bounds/cap');
w.radial_info();const info=read();console.error('Grid',info,'correction groups',count);
for(let first=0;first<count;first++){w.radial_batch(first,1);if(first%100===0)console.error('group',first,'seconds',(performance.now()-t)/1000);}
w.radial_refine(128);const len=w.matches(),pack=read(),best=pack.slice(len);
w.fitted();const fitted=read();
console.error('Finished', (performance.now()-t)/1000,'seconds',best);
console.log(JSON.stringify({info,best,seconds:(performance.now()-t)/1000,scores:pack.slice(0,len),fitted}));
