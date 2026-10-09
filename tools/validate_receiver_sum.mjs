import fs from 'node:fs';
import assert from 'node:assert/strict';
import {prepareReceiverTrains} from '../web/lab/receiver-trains.mjs';
import {projectReceivers} from '../web/lab/beamforming.mjs';
const w=(await WebAssembly.instantiate(fs.readFileSync('web/lab/core.wasm'),{})).instance.exports;
const samples=32,rows=200,pulses=4,fsamp=12.5e6,period=25.37e-6;
const read=()=>new Float32Array(w.memory.buffer,w.result_ptr(),w.result_len()).slice();
let seed=421;const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296-.5;};
const receivers=Array.from({length:4},(_,a)=>Float32Array.from({length:2*samples*rows},(_,i)=>100*(a+1)+(a+1)*random()));
let ptr=w.allocate(receivers[0].length);new Float32Array(w.memory.buffer,ptr,receivers[0].length).set(receivers[0]);
w.load(ptr,rows,samples,100,fsamp,2e-6,period,77e9,99e12);w.release(ptr,receivers[0].length);
w.radial_template(.5,-300,1e5,pulses,0);const q=read();
// Equal target amplitudes with opposite receiver phases cancel in voltage sum.
for(let a=0;a<4;a++)for(let j=0;j<q.length;j++)receivers[a][2*150*samples+j]+=(a%2?-1:1)*20*q[j];
receivers[1].fill(0,2*3*samples,2*4*samples); // Independent mean validity.
const prepared=prepareReceiverTrains(receivers,{samples,rows,perFrame:100,start:150,pulses,bgStart:0,bgStop:20,noiseStart:20,noiseStop:120});
assert.deepEqual(prepared.meanCounts,[20,19,20,20]);
assert(prepared.referenceStarts.every(k=>k>=20&&k+pulses<=120&&Math.floor(k/100)===Math.floor((k+pulses-1)/100)));
for(let a=0;a<4;a++) {
  const first=prepared.referenceStarts[0],offset=(4+a*prepared.quietCount)*q.length;
  assert.deepEqual([...prepared.data.subarray(offset,offset+q.length)],[...receivers[a].subarray(2*first*samples,2*(first+pulses)*samples)],'Noise trains must remain raw I/Q');
}
const beam=projectReceivers(receivers,q,{samples,rows,perFrame:100,start:150,pulses,bgStart:0,bgStop:20,noiseStart:20,noiseStop:120,conjugated:false});
for(let a=0;a<4;a++) {
  const start=prepared.referenceStarts[0];let re=0,im=0;
  for(let j=0;j<q.length;j+=2){const r=receivers[a][2*start*samples+j],i=receivers[a][2*start*samples+j+1];re+=q[j]*r+q[j+1]*i;im+=q[j]*i-q[j+1]*r;}
  assert.equal(beam.noise[0][2*a],re);assert.equal(beam.noise[0][2*a+1],im);
}
w.background(0,20);assert(w.prepare(150,pulses,20,120,0)>=8);
ptr=w.allocate(prepared.data.length);new Float32Array(w.memory.buffer,ptr,prepared.data.length).set(prepared.data);
assert.equal(w.search_receivers(ptr,prepared.data.length,4,prepared.quietCount,2),prepared.quietCount);w.release(ptr,prepared.data.length);
assert.equal(w.search_channels(),4);
const project=train=>{let r=0,i=0;for(let j=0;j<q.length;j+=2){r+=q[j]*train[j]+q[j+1]*train[j+1];i+=q[j]*train[j+1]-q[j+1]*train[j];}return [r,i];};
const energy=train=>project(train).reduce((s,z)=>s+z*z,0),width=q.length;
const events=Array.from({length:4},(_,a)=>prepared.data.subarray(a*width,(a+1)*width));
const observed=events.reduce((s,z)=>s+energy(z),0);
const quiet=Array.from({length:4*prepared.quietCount},(_,i)=>prepared.data.subarray((4+i)*width,(5+i)*width)).reduce((s,z)=>s+energy(z),0)/prepared.quietCount;
assert(w.radial_begin(.5,.5,-300,-300,1e5,1e5,.05,100)>0);w.radial_batch(0,1);w.radial_refine(128);
const length=w.matches(),best=read().slice(length);
assert(Math.abs(best[3]/(observed/quiet)-1)<2e-6,'Exact four-channel numerator and quiet denominator');
const sum=new Float64Array(width);for(const event of events)for(let i=0;i<width;i++)sum[i]+=event[i];
assert(energy(sum)<observed*.001,'Target voltages cancel; incoherent energy remains');
w.train_values();assert.deepEqual([...read()],[...events[2]],'Display RX is independent of shared search objective');
console.log('PASS: four-channel power sum, independent backgrounds, unequal noise powers, common frame-safe quiet trains, cancelling receiver phases and display selection');
