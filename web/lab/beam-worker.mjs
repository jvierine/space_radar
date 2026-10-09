import {phaseSearch, projectReceivers} from './beamforming.mjs?v=20261009gpu2';
const cache=new Map();
async function channel(spec){
 if(cache.has(spec.sha256))return cache.get(spec.sha256);
 const response=await fetch('datasets/test63/'+spec.file);
 if(!response.ok)throw Error('Cannot load all four receiver streams.');
 const bytes=await response.arrayBuffer();
 if(bytes.byteLength!==spec.bytes)throw Error('Receiver stream length mismatch.');
 const digest=[...new Uint8Array(await crypto.subtle.digest('SHA-256',bytes))].map(v=>v.toString(16).padStart(2,'0')).join('');
 if(digest!==spec.sha256)throw Error('Receiver stream SHA-256 mismatch.');
 const data=new Float32Array(bytes);cache.set(spec.sha256,data);return data;
}
let engine;
onmessage=async({data:job})=>{
 try{
  const {meta,match}=job;
  if(meta.transport.length!==4)throw Error('This phase search requires four synchronized receivers.');
  postMessage({id:job.id,type:'progress',text:'Loading and verifying four synchronized receiver streams…'});
  const receivers=await Promise.all(meta.transport.map(channel));
  engine??=(await WebAssembly.instantiateStreaming(fetch('core.wasm?v=20261009gpu2'),{})).instance.exports;
  const p=meta.parameters;
  const ptr=engine.allocate(2*meta.samples);
  new Float32Array(engine.memory.buffer,ptr,2*meta.samples).fill(0);
  engine.load(ptr,1,meta.samples,meta.chirps_per_frame,p.fs,p.T_adc,match.period,p.f_start,p.freq_slope);
  engine.release(ptr,2*meta.samples);
  (match.model==='radial-quadratic' ? engine.radial_template : engine.template_values)(match.best[0],match.best[1],match.best[2],match.pulses,+match.receiver);
  const q=new Float32Array(engine.memory.buffer,engine.result_ptr(),engine.result_len()).slice();
  const projections=projectReceivers(receivers,q,{samples:meta.samples,rows:meta.total_chirps,perFrame:meta.chirps_per_frame,start:match.start,pulses:match.pulses,bgStart:match.bgRange[0],bgStop:match.bgRange[1],noiseStart:match.noiseRange[0],noiseStop:match.noiseRange[1],conjugated:match.best[4]>0});
  const result=phaseSearch(projections.event,projections.noise,job.steps);
  let amplitudeSum=0,energy=0;
  for(let i=0;i<q.length;i+=2){const amplitude=Math.hypot(q[i],q[i+1]);amplitudeSum+=amplitude;energy+=amplitude*amplitude;}
  const effectiveTime=amplitudeSum**2/(p.fs*energy);
  postMessage({id:job.id,type:'result',result:{...result,effectiveTime,referenceStarts:projections.referenceStarts,meanCount:projections.meanCount}});
 }catch(error){postMessage({id:job.id,type:'error',text:error.message});}
};
