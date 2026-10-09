import {GpuRadial,Cancelled} from './gpu-radial.mjs?v=20261009export8';
let cachedGpu=null,disabledReason=null;
const copy=w=>new Float32Array(w.memory.buffer,w.result_ptr(),w.result_len()).slice();
export function importGpuScores(w,scores){
 const ptr=w.allocate(scores.length);
 try{new Float32Array(w.memory.buffer,ptr,scores.length).set(scores);if(w.radial_gpu_scores(ptr,scores.length)<0)throw Error('GPU bank rejected by Rust');}
 finally{w.release(ptr,scores.length);}
}
export function checkGpuScores(w,info,scores,cancelled=()=>false){
 const groups=info[5];let checked=0;
 for(const index of new Set([0,Math.floor(groups/2),groups-1])){
  if(cancelled())throw new Cancelled();
  w.radial_batch(index,1);
  w.radial_gpu_group(index,0);const plan=copy(w),offset=7+2*plan[0]*plan[1]+plan[6];
  const pairs=new Uint32Array(plan.buffer,offset*4,plan[5]*2).slice();
  w.matches();const cpu=new Float32Array(w.memory.buffer,w.result_ptr(),w.result_len());
  for(let j=0;j<pairs.length;j+=2){const i=pairs[j],expected=cpu[i],actual=scores[i];
   if(Number.isFinite(expected)!==Number.isFinite(actual) || Number.isFinite(expected) && Math.abs(actual-expected)>2e-5+5e-4*Math.abs(expected))throw Error('GPU/CPU score check failed');
   if(Number.isFinite(expected))checked++;
  }
 }
 if(!checked)throw Error('No valid nodes in GPU/CPU check');
 return checked;
}
// A GPU failure discards the entire bank and restarts CPU from a clean state.
export async function runRadial({w,info,begin,backend='auto',cancelled=()=>false,progress=()=>{},createGpu=()=>GpuRadial.create(),yieldUI=()=>new Promise(r=>setTimeout(r,0))}){
 let reason=backend==='cpu'?null:disabledReason,gpuSeconds=null;
 if(backend!=='cpu' && info[5]>=4 && !disabledReason){
  try{
   cachedGpu??=await createGpu();if(cancelled())throw new Cancelled();
   const started=performance.now();
   const scores=await cachedGpu.search(w,info,{cancelled,yieldUI,progress:f=>progress(f,'WebGPU')});
   gpuSeconds=(performance.now()-started)/1000;
   progress(1,'Checking WebGPU scores against CPU');await yieldUI();
   checkGpuScores(w,info,scores,cancelled);if(cancelled())throw new Cancelled();
   importGpuScores(w,scores);
   return {backend:'WebGPU',fallbackReason:null,gpuSeconds};
  }catch(error){
   if(error instanceof Cancelled || cancelled())throw new Cancelled();
   reason=error.message;disabledReason=reason;cachedGpu?.destroy();cachedGpu=null;
   progress(0,'CPU fallback: '+reason);await yieldUI();if(cancelled())throw new Cancelled();begin();
  }
 }
 for(let i=0;i<info[5];i++){
  if(cancelled())throw new Cancelled();w.radial_batch(i,1);progress((i+1)/info[5],'CPU (Rust/Wasm)');await yieldUI();
 }
 return {backend:'CPU (Rust/Wasm)',fallbackReason:reason,gpuSeconds};
}
// Test isolation; also useful when a caller explicitly retries after a device reset.
export function resetGpuBackend(){cachedGpu?.destroy();cachedGpu=null;disabledReason=null;}
