import fs from 'node:fs';
import assert from 'node:assert/strict';
import {Worker} from 'node:worker_threads';
const meta=JSON.parse(fs.readFileSync('web/lab/datasets/test63/metadata.json'));
const worker=new Worker(new URL('./gpu_worker_host.mjs',import.meta.url),{workerData:{gpu:false}});
const queue=[],waiters=[];
worker.on('message',m=>{queue.push(m);waiters.splice(0).forEach(f=>f());});
const next=async type=>{for(;;){const i=queue.findIndex(m=>m.type===type||m.type==='error');if(i>=0){const m=queue.splice(i,1)[0];if(m.type==='error')throw Error(m.message);return m;}await new Promise(r=>waiters.push(r));}};
try{
 worker.postMessage({type:'init',meta,rx:0,bgStart:0,bgStop:125});await next('background');
 worker.postMessage({type:'view',start:750,stop:754,component:0,fftSub:true});const sub=await next('view');
 worker.postMessage({type:'view',start:750,stop:754,component:0,fftSub:false});const raw=await next('view');
 const w=(await WebAssembly.instantiate(fs.readFileSync('web/lab/core.wasm'),{})).instance.exports,p=meta.parameters;
 for(const [message,mode] of [[sub,5],[raw,4]]){
  const expected=new Float32Array(4*256);
  for(const spec of [meta.transport[0]]){
   const bytes=fs.readFileSync('web/lab/datasets/test63/'+spec.file),a=new Float32Array(bytes.buffer,bytes.byteOffset,bytes.byteLength/4),ptr=w.allocate(a.length);
   new Float32Array(w.memory.buffer,ptr,a.length).set(a);w.load(ptr,meta.total_chirps,meta.samples,meta.chirps_per_frame,p.fs,p.T_adc,25.37e-6,p.f_start,p.freq_slope);w.release(ptr,a.length);w.background(0,125);w.image(750,754,mode);
   const values=new Float32Array(w.memory.buffer,w.result_ptr(),w.result_len());for(let i=0;i<values.length;i++)expected[i]=values[i];
  }
  assert.deepEqual(message.images[2].values,expected);
 }
 worker.postMessage({type:'rx',rx:2,bgStart:0,bgStop:125});await next('background');
 worker.postMessage({type:'view',start:750,stop:754,component:0,fftSub:true});const other=await next('view');
 assert.notDeepEqual(other.images[2].values,sub.images[2].values);assert.notDeepEqual(other.images[0].values,sub.images[0].values);
 console.log('PASS: raw and mean-subtracted spectra equal the selected receiver FFT; selector changes voltage and spectrum');
}finally{await worker.terminate();}
