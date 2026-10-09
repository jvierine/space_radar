import assert from 'node:assert/strict';
import {phaseSearch,projectReceivers} from '../web/lab/beamforming.mjs';
const event=Float64Array.from([0,72,144,216].flatMap(deg=>[2*Math.cos(deg*Math.PI/180),2*Math.sin(deg*Math.PI/180)]));
const independent=Array.from({length:8},(_,i)=>{const z=new Float64Array(8);z[2*(i%4)]=2;return z;});
const r=phaseSearch(event,independent,10);
assert.deepEqual(r.indices,[8,6,4]);assert.ok(Math.abs(r.peak-16)<1e-10);
assert.ok(r.single.every(v=>Math.abs(v-4)<1e-10));
const same=Float64Array.from([2,0,2,0,2,0,2,0]);
const correlated=Array.from({length:8},()=>Float64Array.from([1,0,1,0,1,0,1,0]));
assert.ok(Math.abs(phaseSearch(same,correlated,10).peak-4)<1e-10,'Correlated receiver noise must not create a fictitious 4× gain');
const samples=10,rows=30;
const receivers=Array.from({length:4},(_,a)=>Float32Array.from({length:rows*samples*2},(_,i)=>i%2?0:100+a));
// Nonzero quiet means; test event subtraction before coherent projection.
for(let a=0;a<4;a++)for(let j=0;j<samples;j++)receivers[a][2*(25*samples+j)]+=a+1;
const q=Float32Array.from({length:2*samples},(_,i)=>i%2?0:1);
let p=projectReceivers(receivers,q,{samples,rows,perFrame:30,start:25,pulses:1,bgStart:0,bgStop:2,noiseStart:2,noiseStop:24,conjugated:false});
assert.deepEqual([...p.event],[10,0,20,0,30,0,40,0]);assert.equal(p.meanCount,2);
// Missing packet on only RX3 invalidates the common selected train.
receivers[3].fill(0,2*25*samples,2*26*samples);
assert.throws(()=>projectReceivers(receivers,q,{samples,rows,perFrame:30,start:25,pulses:1,bgStart:0,bgStop:2,noiseStart:2,noiseStop:24}),/not intact/);
console.log('PASS: known phase recovery, 6.02 dB independent-noise gain, correlated-noise normalization, per-RX background means, common padding rejection');
