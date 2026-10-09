// WebGPU FFT bank. Rust supplies f64 phase corrections and exact bin/index plans.
// Intermediate spectra and the complete ratio bank stay on the GPU until readback.
// Each transform stays in workgroup memory; no global-memory butterfly passes.
function localFft(length) {
 const threads=Math.min(256,length),slots=Math.max(1,length/(2*threads));
 return `var<workgroup> scratch:array<vec2f,${length}>;
fn transform(lane:u32){
 var upper:array<vec2f,${slots}>;var lower:array<vec2f,${slots}>;
 for(var half=1u;half<${length}u;half*=2u){
  for(var t=lane;t<${length/2}u;t+=${threads}u){
   let j=t%half;let base=t/half*(2u*half);let a=scratch[base+j];let b=scratch[base+j+half];
   let angle=-3.141592653589793*f32(j)/f32(half);let w=vec2f(cos(angle),sin(angle));
   let v=vec2f(b.x*w.x-b.y*w.y,b.x*w.y+b.y*w.x);let slot=t/${threads}u;
   upper[slot]=a+v;lower[slot]=a-v;
  }
  workgroupBarrier();
  for(var t=lane;t<${length/2}u;t+=${threads}u){
   let j=t%half;let base=t/half*(2u*half);let slot=t/${threads}u;
   scratch[base+j]=upper[slot];scratch[base+j+half]=lower[slot];
  }
  workgroupBarrier();
 }
}`;
}
const common=`struct Parameters { columns:u32, nodes:u32, orientation:u32, unused:u32 }
@group(0) @binding(0) var<uniform> p:Parameters;`;
export class Cancelled extends Error { constructor(){super('Search cancelled');this.name='Cancelled';} }
export async function bounded(promise,ms,label){
 let timer;try{return await Promise.race([promise,new Promise((_,reject)=>{timer=setTimeout(()=>reject(Error(label+' timed out')),ms);})]);}finally{clearTimeout(timer);}
}
export class GpuRadial {
 static async create(gpu=globalThis.navigator?.gpu) {
  if(!gpu)throw Error('WebGPU is unavailable');
  const adapter=await bounded(gpu.requestAdapter({powerPreference:'high-performance'}),10000,'GPU adapter');
  if(!adapter)throw Error('No WebGPU adapter');
  const device=await bounded(adapter.requestDevice(),10000,'GPU device');
  const engine=new GpuRadial(device);
  engine.adapterInfo={vendor:adapter.info?.vendor,architecture:adapter.info?.architecture,device:adapter.info?.device,description:adapter.info?.description};
  return engine;
 }
 constructor(device){
  this.device=device;this.failure=null;this.cache=new Map();
  device.lost.then(info=>{this.failure=Error('GPU device lost: '+info.message);});
  device.addEventListener('uncapturederror',e=>{this.failure=Error(e.error.message);});
 }
 check(cancelled){if(cancelled())throw new Cancelled();if(this.failure)throw this.failure;}
 async pipeline(code){
  const device=this.device,module=device.createShaderModule({code});
  const info=await module.getCompilationInfo();const errors=info.messages.filter(m=>m.type==='error');
  if(errors.length)throw Error(errors.map(m=>m.message).join('; '));
  return device.createComputePipelineAsync({layout:'auto',compute:{module,entryPoint:'main'}});
 }
 async kernels(n,pulses,nf,ns,trains,channels){
  const key=[n,pulses,nf,ns,trains,channels].join('/');if(this.cache.has(key))return this.cache.get(key);
  const fast=common+localFft(nf)+`
@group(0) @binding(1) var<storage,read> input:array<vec2f>;
@group(0) @binding(2) var<storage,read> correction:array<vec2f>;
@group(0) @binding(3) var<storage,read_write> output:array<vec2f>;
@compute @workgroup_size(${Math.min(256,nf)}) fn main(@builtin(workgroup_id) group:vec3u,@builtin(local_invocation_index) lane:u32){
 let row=group.x;
 for(var local=lane;local<${nf}u;local+=${Math.min(256,nf)}u){
  let j=reverseBits(local)>>${32-Math.log2(nf)}u;var z=vec2f(0.0);
  if(j<${n}u){
   let v=input[row*${n}u+j];let c=correction[(row%${pulses}u)*${n}u+j];
   let qi=select(-c.y,c.y,p.orientation==1u);
   z=vec2f(v.x*c.x-v.y*qi,v.x*qi+v.y*c.x);
  }
  scratch[local]=z;
 }
 workgroupBarrier();transform(lane);
 for(var j=lane;j<${nf}u;j+=${Math.min(256,nf)}u){output[row*${nf}u+j]=scratch[j];}
}`;
  const gather=common+localFft(ns)+`
@group(0) @binding(1) var<storage,read> input:array<vec2f>;
@group(0) @binding(2) var<storage,read> columns:array<u32>;
@group(0) @binding(3) var<storage,read_write> output:array<vec2f>;
@compute @workgroup_size(${Math.min(256,ns)}) fn main(@builtin(workgroup_id) group:vec3u,@builtin(local_invocation_index) lane:u32){
 let col=group.x%p.columns;let train=group.x/p.columns;
 for(var local=lane;local<${ns}u;local+=${Math.min(256,ns)}u){
  let k=reverseBits(local)>>${32-Math.log2(ns)}u;var z=vec2f(0.0);
  if(k<${pulses}u){z=input[(train*${pulses}u+k)*${nf}u+columns[col]];}
  scratch[local]=z;
 }
 workgroupBarrier();transform(lane);
 for(var k=lane;k<${ns}u;k+=${Math.min(256,ns)}u){output[group.x*${ns}u+k]=scratch[k];}
}`;
  const score=common+`
@group(0) @binding(1) var<storage,read> spectra:array<vec2f>;
@group(0) @binding(2) var<storage,read> lookup:array<vec2u>;
@group(0) @binding(3) var<storage,read_write> scores:array<f32>;
@compute @workgroup_size(128) fn main(@builtin(global_invocation_id) id:vec3u){
 let i=id.x;if(i>=p.nodes){return;}
 let node=lookup[i];var event=0.0;var quiet=0.0;
 for(var t=0u;t<${channels}u;t++){let z=spectra[t*p.columns*${ns}u+node.y];event+=dot(z,z);}
 ${trains===channels?"quiet=bitcast<f32>(p.unused);":`for(var t=${channels}u;t<${trains}u;t++){let z=spectra[t*p.columns*${ns}u+node.y];quiet+=dot(z,z)/${(trains-channels)/channels}.0;}`}
 if(quiet>1e-24){let ratio=event/quiet;
  scores[node.x]=max(scores[node.x],ratio);
 }
}`;
  const values=await Promise.all([this.pipeline(fast),this.pipeline(gather),this.pipeline(score)]);
  this.cache.set(key,values);return values;
 }
 async search(wasm,info,{cancelled=()=>false,progress=()=>{},yieldUI=()=>new Promise(r=>setTimeout(r,0))}={}){
  this.check(cancelled);
  const device=this.device,[nr,na,nv,nf,ns,groups]=info,total=nr*na*nv;
  wasm.radial_gpu_trains();const data=new Float32Array(wasm.memory.buffer,wasm.result_ptr(),wasm.result_len()).slice();
  wasm.radial_gpu_group(0,0);let first=new Float32Array(wasm.memory.buffer,wasm.result_ptr(),wasm.result_len()).slice();
  const [n,pulses,,,quiet]=first,channels=wasm.search_channels?.()??1,trains=quiet+channels;
  if(pulses<2)throw Error('GPU backend requires a multi-chirp train');
  const fastCount=trains*pulses*nf,slowCapacity=Math.min(trains*nf*ns,Math.floor(device.limits.maxStorageBufferBindingSize/8));
  const sizes=[data.byteLength,2*n*pulses*4,nf*4,nr*nv*8,fastCount*8,slowCapacity*8,total*4];
  const maximum=Math.max(...sizes);
  if(maximum>device.limits.maxStorageBufferBindingSize || maximum>device.limits.maxBufferSize)throw Error('Search exceeds GPU buffer limits');
  if(trains*pulses>device.limits.maxComputeWorkgroupsPerDimension || Math.ceil(nr*nv/128)>device.limits.maxComputeWorkgroupsPerDimension)throw Error('Search exceeds GPU dispatch limits');
  if(nf*8>device.limits.maxComputeWorkgroupStorageSize)throw Error('FFT exceeds GPU workgroup memory');
  device.pushErrorScope('validation');device.pushErrorScope('out-of-memory');
  const resources=[];
  const buffer=(size,usage)=>{const b=device.createBuffer({size:Math.max(16,size),usage});resources.push(b);return b;};
  const storage=size=>buffer(size,GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST|GPUBufferUsage.COPY_SRC);
  let validationScope=true,oomScope=true;
  try{
   const [prepare,gather,score]=await bounded(this.kernels(n,pulses,nf,ns,trains,channels),30000,'GPU shader compilation');
   const input=storage(data.byteLength),correction=storage(2*n*pulses*4),columns=storage(nf*4),lookup=storage(nr*nv*8);
   const fastBuffer=storage(fastCount*8),slowBuffer=storage(slowCapacity*8),scores=storage(total*4);
   const readback=buffer(total*4,GPUBufferUsage.MAP_READ|GPUBufferUsage.COPY_DST);
   device.queue.writeBuffer(input,0,data);
   device.queue.writeBuffer(scores,0,new Float32Array(total).fill(-1));
   const group=(pipeline,buffers)=>device.createBindGroup({layout:pipeline.getBindGroupLayout(0),entries:buffers.map((b,binding)=>({binding,resource:{buffer:b}}))});
   const params=(a)=>{const b=buffer(16,GPUBufferUsage.UNIFORM|GPUBufferUsage.COPY_DST);device.queue.writeBuffer(b,0,new Uint32Array(a));return b;};
   const mainParams=params([0,0,0,0]);
   const prepGroup=group(prepare,[mainParams,input,correction,fastBuffer]);
   const gatherGroup=group(gather,[mainParams,fastBuffer,columns,slowBuffer]);
   const scoreGroup=group(score,[mainParams,slowBuffer,lookup,scores]);
   const run=(encoder,pipeline,bind,groups)=>{
    if(!groups)return;
    const pass=encoder.beginComputePass();pass.setPipeline(pipeline);pass.setBindGroup(0,bind);pass.dispatchWorkgroups(groups);pass.end();
   };
   for(let index=0;index<groups;index++){
    this.check(cancelled);
    for(let orientation=0;orientation<2;orientation++){
     wasm.radial_gpu_group(index,orientation);const pack=new Float32Array(wasm.memory.buffer,wasm.result_ptr(),wasm.result_len()).slice();
     const nodes=pack[5],nc=pack[6];if(!nodes)continue;
     if(trains*nc*ns>slowCapacity || trains*nc>device.limits.maxComputeWorkgroupsPerDimension)throw Error("Search exceeds GPU compact-spectrum limits");
     let offset=7;
     device.queue.writeBuffer(correction,0,pack.subarray(offset,offset+=2*n*pulses));
     device.queue.writeBuffer(columns,0,pack.subarray(offset,offset+=nc));
     device.queue.writeBuffer(lookup,0,pack.subarray(offset));
     const parameters=new Uint32Array([nc,nodes,orientation,0]);new Float32Array(parameters.buffer)[3]=wasm.radial_noise_power?.()??0;
     device.queue.writeBuffer(mainParams,0,parameters);
     const encoder=device.createCommandEncoder();run(encoder,prepare,prepGroup,trains*pulses);
     run(encoder,gather,gatherGroup,trains*nc);
     run(encoder,score,scoreGroup,Math.ceil(nodes/128));device.queue.submit([encoder.finish()]);
    }
    // Bound outstanding work so cancellation and device errors are seen promptly.
    if(index%4===3 || index===groups-1){await bounded(device.queue.onSubmittedWorkDone(),30000,'GPU FFT batch');this.check(cancelled);progress((index+1)/groups);await yieldUI();this.check(cancelled);}
   }
   const encoder=device.createCommandEncoder();encoder.copyBufferToBuffer(scores,0,readback,0,total*4);device.queue.submit([encoder.finish()]);
   await bounded(readback.mapAsync(GPUMapMode.READ),30000,'GPU score readback');this.check(cancelled);
   const result=new Float32Array(readback.getMappedRange(),0,total).slice();readback.unmap();
   for(let i=0;i<result.length;i++)if(result[i]===-1)result[i]=NaN;
   const oom=await device.popErrorScope();oomScope=false;
   const validation=await device.popErrorScope();validationScope=false;
   if(oom || validation)throw Error((oom || validation).message);
   if(!result.some(Number.isFinite) || result.some(v=>v===Infinity || v===-Infinity || v<0))throw Error('GPU returned invalid matched-energy scores');
   return result;
  }finally{
   for(const b of resources)b.destroy();
   if(oomScope)await device.popErrorScope();if(validationScope)await device.popErrorScope();
  }
 }
 destroy(){this.device.destroy();}
}
