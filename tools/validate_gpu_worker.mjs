import fs from 'node:fs';
import assert from 'node:assert/strict';
import {Worker} from 'node:worker_threads';
const meta=JSON.parse(fs.readFileSync('web/lab/datasets/test63/metadata.json'));
const job={type:'search',algorithm:'fft',computeBackend:'auto',start:750,pulses:8,period:25.37e-6,receiver:false,noiseStart:125,noiseStop:625,grid:{xMin:.28,xMax:.32,yMin:3e5,yMax:7e5,vMin:280,vMax:360},loss:5};
async function session(gpu){
 const worker=new Worker(new URL('./gpu_worker_host.mjs',import.meta.url),{workerData:{gpu}});
 const queue=[],waiters=[];let failure;
 worker.on('message',m=>{queue.push(m);waiters.splice(0).forEach(w=>w());});
 worker.on('error',e=>{failure=e;waiters.splice(0).forEach(w=>w());});
 const next=async type=>{
  const timeout=setTimeout(()=>{failure=Error('Worker test timeout');waiters.splice(0).forEach(w=>w());},60000);
  try{for(;;){if(failure)throw failure;const i=queue.findIndex(m=>m.type==='error'||m.type===type);if(i>=0){const m=queue.splice(i,1)[0];if(m.type==='error')throw Error(m.message);return m;}await new Promise(r=>waiters.push(r));}}finally{clearTimeout(timeout);}
 };
 worker.postMessage({type:'init',meta,rx:0,bgStart:0,bgStop:125});await next('background');
 worker.postMessage({type:'period',period:job.period});
 return {worker,next,queue};
}
let cpu;
for(const gpu of [false,true]){
 const s=await session(gpu);
 try{
  s.worker.postMessage(job);const result=(await s.next('match')).result;await s.next('complete');
  assert.equal(result.backend,gpu?'WebGPU':'CPU (Rust/Wasm)');
  if(!gpu){assert.match(result.fallbackReason,/unavailable/);cpu=result;}
  else {assert.deepEqual(result.best,cpu.best);assert.equal(result.cube.length,cpu.cube.length);}
  s.queue.length=0;
  // Cancellation followed immediately by another search must not mix banks.
  s.worker.postMessage({...job,grid:{xMin:.001,xMax:3,yMin:0,yMax:1e6,vMin:-900,vMax:900}});
  const progress=await s.next('progress');assert(progress.fraction>=0);
  s.worker.postMessage({type:'cancel'});s.worker.postMessage({...job,start:875});await s.next('cancelled');
  const replacement=(await s.next('match')).result;await s.next('complete');assert.equal(replacement.start,875);assert.equal(replacement.cube.length,cpu.cube.length);
  assert(!s.queue.some(m=>m.type==='match' && m.result.start===750),'Cancelled bank must not publish a fit');
  console.log('PASS: browser worker protocol, '+(gpu?'WebGPU':'unavailable-GPU CPU fallback')+', cancellation and immediate replacement search');
 }finally{await s.worker.terminate();}
}
