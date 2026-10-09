import {prepareReceiverTrains} from '../web/lab/receiver-trains.mjs';
// Native Dawn/Metal executes the same WebGPU shaders used by the worker.
// npm install --prefix /tmp/fmcw-webgpu-test webgpu
// WEBGPU_MODULE=/tmp/fmcw-webgpu-test/node_modules/webgpu/index.js node tools/validate_gpu_radial.mjs
import fs from 'node:fs';
import assert from 'node:assert/strict';
import {pathToFileURL} from 'node:url';
import {GpuRadial,Cancelled,bounded} from '../web/lab/gpu-radial.mjs?v=20261009drag3';
import {runRadial,resetGpuBackend,checkGpuScores,importGpuScores} from '../web/lab/radial-backend.mjs';
const {create,globals}=await import(pathToFileURL(process.env.WEBGPU_MODULE??'/tmp/fmcw-webgpu-test/node_modules/webgpu/index.js'));
Object.assign(globalThis,globals);
const provider=create(['backend=metal']),gpu=await GpuRadial.create(provider);
const read=w=>new Float32Array(w.memory.buffer,w.result_ptr(),w.result_len()).slice();
const p={rows:1000,n:225,frame:125,fs:12.5e6,adc:2e-6,period:25.37e-6,f0:77e9,slope:99.98738765716553e12};
function measurement(pulses,truth,orientation){
 let seed=7539;const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return(seed+.5)/4294967296;};
 const data=new Float32Array(p.rows*p.n*2);
 for(let i=0;i<data.length;i+=2){const radius=20*Math.sqrt(-2*Math.log(random())),angle=2*Math.PI*random();data[i]=radius*Math.cos(angle)+400;data[i+1]=radius*Math.sin(angle)-100;}
 const centre=p.adc+(p.n-1)/(2*p.fs)+(pulses-1)*p.period/2;
 for(let k=0;k<pulses;k++)for(let j=0;j<p.n;j++){
  const u=p.adc+j/p.fs,h=k*p.period+u-centre,r=truth[0]+truth[1]*h+.5*truth[2]*h*h,tau=2*r/299792458;
  const phase=2*Math.PI*(-(p.f0+p.slope*u)*tau+.5*p.slope*tau*tau),i=2*((750+k)*p.n+j);
  data[i]+=100*Math.cos(phase);data[i+1]+=(orientation?-1:1)*100*Math.sin(phase);
 }
 // Padding controls must be excluded by Rust before either backend sees trains.
 data.fill(0,2*880*p.n,2*884*p.n);return data;
}
async function fixture(pulses,bounds,truth,orientation=0,hpf=0,channels=1){
 const w=(await WebAssembly.instantiate(fs.readFileSync('web/lab/core.wasm'),{})).instance.exports;
 const bytes=truth?null:fs.readFileSync('web/lab/datasets/test63/rx0.f32');
 const data=truth?measurement(pulses,truth,orientation):new Float32Array(bytes.buffer,bytes.byteOffset,bytes.byteLength/4),ptr=w.allocate(data.length);new Float32Array(w.memory.buffer,ptr,data.length).set(data);
 w.load(ptr,data.length/(2*p.n),p.n,p.frame,p.fs,p.adc,p.period,p.f0,p.slope);w.release(ptr,data.length);w.background(0,125);assert(w.prepare(750,pulses,125,625,hpf)>=8);
 if(channels===4){
  const receivers=Array.from({length:4},(_,i)=>{
   if(!truth){const b=fs.readFileSync(`web/lab/datasets/test63/rx${i}.f32`);return new Float32Array(b.buffer,b.byteOffset,b.byteLength/4);}
   const d=measurement(pulses,truth,orientation),angle=i*.7,scale=1+i*.15;
   for(let j=0;j<d.length;j+=2){const r=d[j],im=d[j+1];d[j]=scale*(Math.cos(angle)*r-Math.sin(angle)*im);d[j+1]=scale*(Math.sin(angle)*r+Math.cos(angle)*im);}
   return d;
  });
  const prepared=prepareReceiverTrains(receivers,{samples:p.n,rows:receivers[0].length/(2*p.n),perFrame:p.frame,start:750,pulses,bgStart:0,bgStop:125,noiseStart:125,noiseStop:625});
  const ptr=w.allocate(prepared.data.length);new Float32Array(w.memory.buffer,ptr,prepared.data.length).set(prepared.data);assert(w.search_receivers(ptr,prepared.data.length,4,prepared.quietCount,0)>=8);w.release(ptr,prepared.data.length);
 }
 const begin=()=>{const groups=w.radial_begin(...bounds,.05,32000000);assert(groups>0);return groups;};
 begin();w.radial_info();return {w,begin,info:Array.from(read(w))};
}
function difference(a,b){
 assert.equal(a.length,b.length);let maximum=0;
 for(let i=0;i<a.length;i++){
  assert.equal(Number.isFinite(a[i]),Number.isFinite(b[i]),'Valid-node mask at '+i);
  if(Number.isFinite(a[i])){const error=Math.abs(a[i]-b[i]);maximum=Math.max(maximum,error/(1+Math.abs(b[i])));assert(error<=2e-5+2e-4*Math.abs(b[i]),`Node ${i}: GPU ${a[i]}, CPU ${b[i]}`);}
 }
 return maximum;
}
function projections(scores,[nr,na,nv]){
 const vr=new Float32Array(nr*nv).fill(NaN),va=new Float32Array(na*nv).fill(NaN);
 for(let ir=0;ir<nr;ir++)for(let ia=0;ia<na;ia++)for(let iv=0;iv<nv;iv++){
  const value=scores[(ir*na+ia)*nv+iv];if(!Number.isFinite(value))continue;
  for(const [array,i] of [[vr,ir*nv+iv],[va,ia*nv+iv]])if(!Number.isFinite(array[i]) || value>array[i])array[i]=value;
 }
 return [vr,va];
}
function finish(w){assert(w.radial_refine(128)>=0);const n=w.matches(),a=read(w);return {scores:a.slice(0,n),best:a.slice(n)};}
const records=[];
try{
 for(const [name,pulses,bounds,truth,orientation,hpf,channels] of [
  ['wide_8',8,[.001,3,-900,900,0,1e6],[.3,-320,5e5],0,0],
  ['aliased_conjugate_8',8,[.28,.32,280,360,3e5,7e5],[.3,320,5e5],1,0],
  ['focused_16_hpf',16,[.28,.32,-340,-300,4e5,6e5],[.3,-320,5e5],0,1],
  ['signed_acceleration_4',4,[.02,.06,-500,500,-1e6,1e6],[.04,-320,-5e5],0,0],
  ['singleton_2',2,[.3,.3,-320,-320,5e5,5e5],[.3,-320,5e5],0,0],
  ['test63_8',8,[.001,3,-900,900,0,1e6],null,null,0],
  ['four_rx_8',8,[.28,.32,-340,-300,4e5,6e5],[.3,-320,5e5],0,0,4],
  ['test63_four_rx_8',8,[.001,3,-900,900,0,1e6],null,null,0,4]
 ]){
  const {w,begin,info}=await fixture(pulses,bounds,truth,orientation,hpf,channels);
  const start=performance.now();const scores=await gpu.search(w,info);const cold=(performance.now()-start)/1000;
  const checkStart=performance.now();checkGpuScores(w,info,scores);importGpuScores(w,scores);const gpuResult=finish(w);const checkTime=(performance.now()-checkStart)/1000;
  begin();const cpuStart=performance.now();for(let i=0;i<info[5];i++)w.radial_batch(i,1);const cpu=finish(w),cpuTime=(performance.now()-cpuStart)/1000;
  const maxError=difference(scores,cpu.scores),gp=projections(scores,info),cp=projections(cpu.scores,info);difference(gp[0],cp[0]);difference(gp[1],cp[1]);
  assert.deepEqual([...gpuResult.best.slice(0,3)],[...cpu.best.slice(0,3)]);assert.equal(gpuResult.best[4],cpu.best[4]);if(orientation!==null)assert.equal(cpu.best[4],orientation);assert(Math.abs(gpuResult.best[3]/cpu.best[3]-1)<1e-6);
  begin();const warmStart=performance.now();const warm=await gpu.search(w,info);checkGpuScores(w,info,warm);importGpuScores(w,warm);finish(w);const warmTime=(performance.now()-warmStart)/1000;
  const row={name,pulses,bounds,grid:info,coldSeconds:cold+checkTime,warmSeconds:warmTime,cpuSeconds:cpuTime,speedup:cpuTime/warmTime,maxNormalizedScoreError:maxError,best:[...cpu.best]};records.push(row);
  console.error(JSON.stringify(row));
 }
 // Fault-injection uses the real Rust bank; verify complete clean CPU recovery.
 const {w,begin,info}=await fixture(8,[.28,.32,280,360,3e5,7e5],[.3,320,5e5],1);
 for(let i=0;i<info[5];i++)w.radial_batch(i,1);const expected=finish(w).scores;
 assert.throws(()=>importGpuScores(w,new Float32Array(expected.length+1)),/rejected/);
 const negative=expected.slice();negative[0]=-3;assert.throws(()=>importGpuScores(w,negative),/rejected/);
 w.matches();difference(read(w).slice(0,expected.length),expected);
 for(const reason of ['No WebGPU adapter','GPU device lost','GPU buffer allocation failed','GPU shader compilation failed','GPU/CPU score check failed']){
  resetGpuBackend();begin();
 const badScores=expected.map(v=>Number.isFinite(v)?v*2:v);
 const bad=await runRadial({w,info,begin,createGpu:async()=>({destroy(){},search:async()=>badScores}),yieldUI:async()=>{}});
 assert.equal(bad.backend,'CPU (Rust/Wasm)');assert.match(bad.fallbackReason,/score check/);
 let checks=0;begin();await assert.rejects(gpu.search(w,info,{cancelled:()=>++checks>3}),Cancelled);
 resetGpuBackend();begin();let resets=0;
  const result=await runRadial({w,info,begin:()=>{resets++;begin();},createGpu:async()=>({destroy(){},async search(){w.radial_batch(0,1);throw Error(reason);}}),yieldUI:async()=>{}});
  assert.equal(result.backend,'CPU (Rust/Wasm)');assert.equal(result.fallbackReason,reason);assert.equal(resets,1);w.matches();difference(read(w).slice(0,expected.length),expected);
  begin();await runRadial({w,info,begin,createGpu:()=>{throw Error('Must not retry failed GPU');},yieldUI:async()=>{}});
 }
 resetGpuBackend();begin();
 const badScores=expected.map(v=>Number.isFinite(v)?v*2:v);
 const bad=await runRadial({w,info,begin,createGpu:async()=>({destroy(){},search:async()=>badScores}),yieldUI:async()=>{}});
 assert.equal(bad.backend,'CPU (Rust/Wasm)');assert.match(bad.fallbackReason,/score check/);
 let checks=0;begin();await assert.rejects(gpu.search(w,info,{cancelled:()=>++checks>3}),Cancelled);
 resetGpuBackend();begin();let resets=0;
 await assert.rejects(runRadial({w,info,begin:()=>{resets++;begin();},cancelled:()=>true,createGpu:async()=>({destroy(){}})}),Cancelled);assert.equal(resets,0,'Cancellation must not start CPU');
 resetGpuBackend();begin();const forced=await runRadial({w,info,begin,backend:'cpu',createGpu:()=>{throw Error('Must not request GPU');},yieldUI:async()=>{}});assert.equal(forced.backend,'CPU (Rust/Wasm)');
 // Actual device loss, plus real limits check before allocation.
 const lost=await GpuRadial.create(provider);lost.destroy();await lost.device.lost;await assert.rejects(lost.search(w,info),/lost/);
 await assert.rejects(gpu.search(w,[50000000,1,1,...info.slice(3)]),/buffer limits/);
 await assert.rejects(bounded(new Promise(()=>{}),5,'Test GPU hang'),/timed out/);
 console.error('PASS: full grids, masks, both MAX maps, exact peaks, both I/Q orientations, HPF refinement, cancellation, device loss, limits and CPU recovery');
 console.log(JSON.stringify({generator:'tools/validate_gpu_radial.mjs',gpu:{...gpu.adapterInfo,backend:'Dawn/Metal',platform:process.platform,arch:process.arch},records}));
}finally{resetGpuBackend();gpu.destroy();}
