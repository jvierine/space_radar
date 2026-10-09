import fs from 'node:fs';import assert from 'node:assert/strict';
import {prepareReceiverTrains} from '../web/lab/receiver-trains.mjs';
import {projectReceivers,phaseSearch} from '../web/lab/beamforming.mjs';
const meta=JSON.parse(fs.readFileSync('web/lab/datasets/test63/metadata.json')),p=meta.parameters;
const receivers=meta.transport.map(s=>{const b=fs.readFileSync('web/lab/datasets/test63/'+s.file);return new Float32Array(b.buffer,b.byteOffset,b.byteLength/4);});
const options={samples:225,rows:6250,perFrame:125,start:3755,pulses:8,bgStart:3648,bgStop:3649,noiseStart:3648,noiseStop:3649,fullBandwidth:true};
const prepared=prepareReceiverTrains(receivers,options);assert.equal(prepared.noiseSamples,225);
for(let a=0;a<4;a++){let power=0;for(let j=0;j<225;j++){const i=2*(3648*225+j);power+=receivers[a][i]**2+receivers[a][i+1]**2;}assert(Math.abs(prepared.noisePower[a]/(power/225)-1)<1e-6);}
const w=(await WebAssembly.instantiate(fs.readFileSync('web/lab/core.wasm'),{})).instance.exports,read=()=>new Float32Array(w.memory.buffer,w.result_ptr(),w.result_len()).slice();
let ptr=w.allocate(receivers[0].length);new Float32Array(w.memory.buffer,ptr,receivers[0].length).set(receivers[0]);w.load(ptr,6250,225,125,p.fs,p.T_adc,25.37e-6,p.f_start,p.freq_slope);w.release(ptr,receivers[0].length);w.background(3648,3649);
assert.equal(w.prepare_power(3755,8,3648,3649,0),1);
ptr=w.allocate(prepared.data.length);new Float32Array(w.memory.buffer,ptr,prepared.data.length).set(prepared.data);assert.equal(w.search_receivers(ptr,prepared.data.length,4,1,0),1);w.release(ptr,prepared.data.length);
ptr=w.allocate(4);new Float32Array(w.memory.buffer,ptr,4).set(prepared.noisePower);assert.equal(w.search_noise_power(ptr,4),4);w.release(ptr,4);
w.radial_begin(.5,.5,320,320,5e5,5e5,.05,100);w.radial_batch(0,1);w.radial_refine(128);const n=w.matches(),best=read().slice(n);
w.radial_template(...best.slice(0,3),8,0);const q=read(),energy=q.reduce((s,v)=>s+v*v,0);let signal=0;
for(let a=0;a<4;a++){let re=0,im=0;const z=prepared.data.subarray(a*q.length,(a+1)*q.length);for(let j=0;j<q.length;j+=2){const qi=best[4]?-q[j+1]:q[j+1];re+=q[j]*z[j]+qi*z[j+1];im+=q[j]*z[j+1]-qi*z[j];}signal+=re*re+im*im;}
assert(Math.abs(best[3]/(signal/(prepared.noisePower.reduce((a,b)=>a+b,0)*energy))-1)<1e-5);
const projected=projectReceivers(receivers,q,{...options,conjugated:best[4]>0}),beam=phaseSearch(projected.event,projected.noise,10,projected.noiseCovariance);
let rawBeamPower=0;for(let j=0;j<225;j++){let re=0,im=0;for(let a=0;a<4;a++){const i=2*(3648*225+j),c=Math.cos(beam.phases[a])/2,s=Math.sin(beam.phases[a])/2;re+=c*receivers[a][i]-s*receivers[a][i+1];im+=s*receivers[a][i]+c*receivers[a][i+1];}rawBeamPower+=re*re+im*im;}
assert(Math.abs(beam.reference/(rawBeamPower/225*energy)-1)<1e-8);
console.log('PASS: one raw chirp, full-bandwidth complex sample power, exact matched-filter noise scaling, and beamformed raw covariance normalization');
