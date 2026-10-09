// Regression for the GUI state reported by the user: 29 background chirps, 8-chirp integration.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import {Worker} from 'node:worker_threads';
import {prepareReceiverTrains,commonValid} from '../web/lab/receiver-trains.mjs';
const meta=JSON.parse(fs.readFileSync('web/lab/datasets/test63/metadata.json'));
const receivers=meta.transport.map(spec=>{const b=fs.readFileSync('web/lab/datasets/test63/'+spec.file);return new Float32Array(b.buffer,b.byteOffset,b.byteLength/4);});
const valid=commonValid(receivers,meta.samples,meta.total_chirps),expected=[];
for(let start=3701;start+8<=3916;start++)if(Math.floor(start/125)===Math.floor((start+7)/125)&&valid.slice(start,start+8).every(Boolean))expected.push(start);
assert.equal(expected.length,194);
const opts={samples:225,rows:6250,perFrame:125,start:3701,pulses:8,bgStart:3648,bgStop:3677,noiseStart:3648,noiseStop:3677,fullBandwidth:true};
const prepared=prepareReceiverTrains(receivers,opts);assert.equal(prepared.quietCount,1);assert.equal(prepared.noiseSamples,29*225);
assert.deepEqual(prepared.referenceStarts,Array.from({length:29},(_,i)=>3648+i));
const one=prepareReceiverTrains(receivers,{...opts,bgStop:3649,noiseStop:3649});assert.equal(one.noiseSamples,225);
assert.throws(()=>prepareReceiverTrains(receivers,{...opts,bgStart:3750,bgStop:3751,noiseStart:3750,noiseStop:3751,start:3750}),/one intact chirp/);
const worker=new Worker(new URL('./gpu_worker_host.mjs',import.meta.url),{workerData:{gpu:true}}),queue=[],waiters=[];
let failure;worker.on('message',m=>{queue.push(m);waiters.splice(0).forEach(f=>f());});worker.on('error',e=>{failure=e;waiters.splice(0).forEach(f=>f());});
const next=async type=>{
 const timer=setTimeout(()=>{failure=Error('Short-background test timed out');waiters.splice(0).forEach(f=>f());},60000);
 try{for(;;){if(failure)throw failure;const i=queue.findIndex(m=>m.type===type||m.type==='error');if(i>=0)return queue.splice(i,1)[0];await new Promise(r=>waiters.push(r));}}finally{clearTimeout(timer);}
};
const job={type:'search',algorithm:'fft',computeBackend:'cpu',start:3755,pulses:8,period:25.37e-6,receiver:false,noiseStart:3648,noiseStop:3677,grid:{xMin:.5,xMax:.5,yMin:5e5,yMax:5e5,vMin:320,vMax:320},loss:5,scan:true,scanStart:3701,scanStop:3916,scanStride:1,beamSteps:3};
try{
 worker.postMessage({type:'init',meta,rx:0,bgStart:3648,bgStop:3677});assert.equal((await next('background')).type,'background');
 worker.postMessage({type:'period',period:job.period});worker.postMessage(job);
 const complete=await next('complete');assert.equal(complete.type,'complete',complete.message);
 const starts=queue.filter(m=>m.type==='match'&&m.result.scanPoint).map(m=>m.result.start);assert.deepEqual(starts,expected);
 assert(queue.filter(m=>m.type==='match').every(m=>Number.isFinite(m.result.beam.peak)));
 queue.length=0;
 // The exact full user grid must also complete with those short background bounds.
 worker.postMessage({...job,computeBackend:'auto',grid:{xMin:.001,xMax:3,yMin:0,yMax:1e6,vMin:0,vMax:900},scanStop:3709});
 const full=await next('complete');assert.equal(full.type,'complete',full.message);assert.equal(full.results[0].backend,'WebGPU');console.log('User settings peak: '+(10*Math.log10(Math.max(full.results[0].best[3]-1,1e-12))).toFixed(2)+' dB; background residual RX powers: '+Array.from(full.results[0].noisePower).join(', '));assert.equal(full.results[0].referenceStarts.length,29);assert(full.results[0].best[3]>2,'Shared-view peak has positive excess-signal SNR');assert(Math.abs(full.results[0].beam.analysisBandwidth-12.5e6/1800)<.01);for(let a=0;a<4;a++)assert(Math.abs(full.results[0].beam.analysisNoisePower[a]/(full.results[0].beam.rawNoisePower[a]/1800)-1)<1e-6);
 queue.length=0;
 worker.postMessage({...job,computeBackend:'auto',grid:{xMin:.001,xMax:3,yMin:0,yMax:1e6,vMin:0,vMax:900},scanStart:3755,scanStop:3763});
 const target=await next('complete');assert.equal(target.type,'complete',target.message);
 console.log('Selected chirp 3755 RX diagnostics: '+JSON.stringify(target.results[0].beam.channelPower.map((p,i)=>({rx:i,noise:target.results[0].beam.rawNoisePower[i],analysisNoise:target.results[0].beam.analysisNoisePower[i],snrDb:10*Math.log10(Math.max(target.results[0].beam.single[i]-1,1e-12)),matchedPower:p.observed}))));
 queue.length=0;worker.postMessage({type:'background',start:3648,stop:3649});await next('background');
 worker.postMessage({...job,noiseStop:3649,scanStop:3709});const single=await next('complete');assert.equal(single.type,'complete',single.message);assert.equal(single.results[0].noiseSamples,225);
 assert(single.results[0].beam.single.every(Number.isFinite));
 console.log('PASS: exact user background, all 194 analysis starts, full user WebGPU grid, finite beam outputs, and one-chirp full-bandwidth noise estimate');
}finally{await worker.terminate();}
