import fs from 'node:fs';
import assert from 'node:assert/strict';
import {Worker} from 'node:worker_threads';
const worker=new Worker(new URL('./gpu_worker_host.mjs',import.meta.url),{workerData:{gpu:false}});
const meta=JSON.parse(fs.readFileSync('web/lab/datasets/test63/metadata.json'));
let resolve,reject;let expectedType='background';const points=[];
const timer=setTimeout(()=>reject(Error('Scan test timeout')),60000);
const done=new Promise((r,j)=>{resolve=r;reject=j;});
worker.on('error',reject);
worker.on('message',m=>{
 try {
  if(m.type==='error')throw Error(m.message);
  if(m.type==='background'&&expectedType==='background'){
   expectedType='complete';worker.postMessage({type:'period',period:25.37e-6});
   worker.postMessage({type:'search',algorithm:'fft',computeBackend:'cpu',start:750,pulses:8,period:25.37e-6,receiver:false,noiseStart:0,noiseStop:125,grid:{xMin:.3,xMax:.3,yMin:5e5,yMax:5e5,vMin:-320,vMax:-320},loss:5,scan:true,scanStart:750,scanStop:762,scanStride:10,beamSteps:3});
  }
  if(m.type==='match'&&m.result.scanPoint){
   const r=m.result;points.push(r.start);assert.equal(r.searchReceivers,4);assert.equal(r.beam.single.length,4);assert.equal(r.beam.phases.length,4);assert(Number.isFinite(r.beam.peak));
  }
  if(m.type==='complete'){
   assert.deepEqual(points,[750,751,752,753,754],'Every contained eight-chirp train, even with obsolete stride input');
   assert.equal(m.results[0].scanResults.length,25);resolve();
  }
 }catch(e){reject(e);}
});
worker.postMessage({type:'init',meta,rx:0,bgStart:0,bgStop:125});
try{await done;console.log('PASS: every intact train start, incremental matches and beam/SNR results before scan completion');}finally{clearTimeout(timer);await worker.terminate();}
