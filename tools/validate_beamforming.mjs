import assert from 'node:assert/strict';
import {phaseSearch,projectReceivers,receiverGrowth} from '../web/lab/beamforming.mjs';
const event=Float64Array.from([0,72,144,216].flatMap(deg=>[2*Math.cos(deg*Math.PI/180),2*Math.sin(deg*Math.PI/180)]));
const independent=Array.from({length:8},(_,i)=>{const z=new Float64Array(8);z[2*(i%4)]=2;return z;});
const r=phaseSearch(event,independent,10);
assert.deepEqual(r.indices,[8,6,4]);assert.ok(Math.abs(r.peak-16)<1e-10);
assert.ok(r.single.every(v=>Math.abs(v-4)<1e-10));
assert.ok(Math.abs(r.total.singleAverage-4)<1e-10);
assert.ok(Math.abs(r.total.gain-4)<1e-10);
assert.ok(Math.abs(r.total.percentIdeal-100)<1e-8);
for(const g of r.growth) {
  assert.ok(Math.abs(g.score-4*g.count)<1e-9);
  if(g.count>1)assert.ok(Math.abs(g.percentIdeal-100)<1e-8);
}
// Independently evaluate every phase triple and MAX onto each axis pair.
const expected=Array.from({length:3},()=>new Float32Array(100).fill(-Infinity));
for(let i=0;i<10;i++)for(let j=0;j<10;j++)for(let k=0;k<10;k++) {
  const phases=[0,i,j,k].map(x=>x*2*Math.PI/10);
  const energy=z=>{
    const re=z.reduce((s,_,a)=>a%2?s:s+(z[a]*Math.cos(phases[a/2])-z[a+1]*Math.sin(phases[a/2]))/2,0);
    const im=z.reduce((s,_,a)=>a%2?s:s+(z[a]*Math.sin(phases[a/2])+z[a+1]*Math.cos(phases[a/2]))/2,0);
    return re*re+im*im;
  };
  const value=10*Math.log10(Math.max(energy(event)/(independent.reduce((s,z)=>s+energy(z),0)/independent.length)-1,1));
  for(const [p,index] of [[0,j*10+i],[1,k*10+i],[2,k*10+j]])expected[p][index]=Math.max(expected[p][index],value);
}
assert.deepEqual(r.projections,expected);
for(const [p,[x,y]] of [[0,[0,1]],[1,[0,2]],[2,[1,2]]])
  assert.ok(Math.abs(r.projections[p][r.indices[y]*10+r.indices[x]]-10*Math.log10(Math.max(r.gridPeak-1,1)))<1e-5);
console.log('PASS: all three MAX phase projections, axis ordering, and projected grid-peak cells');
const same=Float64Array.from([2,0,2,0,2,0,2,0]);
const correlated=Array.from({length:8},()=>Float64Array.from([1,0,1,0,1,0,1,0]));
assert.ok(Math.abs(phaseSearch(same,correlated,10).peak-4)<1e-10,'Correlated receiver noise must not create a fictitious 4× gain');
for(const g of phaseSearch(same,correlated,10).growth) {
  assert.ok(Math.abs(g.score-4)<1e-9);
  if(g.count>1)assert.ok(Math.abs(g.percentIdeal)<1e-8);
}
assert.ok(Math.abs(phaseSearch(same,correlated,10).total.percentIdeal-25)<1e-8);
const unequal=phaseSearch([2,0,1,0,.5,0,.25,0],independent,10);
assert.ok(Math.abs(unequal.total.singleAverage-(4+1+.25+.0625)/4)<1e-10);
assert.ok(Math.abs(unequal.total.gain-unequal.peak/unequal.total.singleAverage)<1e-10);
const weak=receiverGrowth([2,0,0,0,0,0,0,0],independent,[0,0,0,0]);
for(const g of weak) {
  assert.ok(Math.abs(g.score-4/g.count)<1e-10);
  if(g.count>1){assert.ok(g.gain<1);assert.ok(!Number.isFinite(g.percentIdeal));}
}
console.log('PASS: receiver increments at 100% ideal, correlated-noise zero gain, and negative contributions');
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

// Off-grid phases and a solution straddling the periodic boundary.
const planted=[0,.11,2.3,2*Math.PI-.023];
const offGrid=Float64Array.from(planted.flatMap(p=>[2*Math.cos(p),2*Math.sin(p)]));
const refined=phaseSearch(offGrid,independent,10);
assert.ok(refined.refinement.converged);
assert.ok(refined.peak>=refined.gridPeak);
assert.ok(refined.peak>refined.gridPeak+.01,'Off-grid peak must improve');
for(let j=1;j<4;j++)assert.ok(Math.abs(Math.atan2(Math.sin(refined.phases[j]+planted[j]),Math.cos(refined.phases[j]+planted[j])))<1e-5);
assert.ok(Math.abs(refined.peak-16)<1e-9);
assert.ok(Math.abs(refined.growth[3].score-refined.peak)<1e-10);
const combined=z=>{
 let re=0,im=0;for(let j=0;j<4;j++){const c=Math.cos(refined.phases[j])/2,s=Math.sin(refined.phases[j])/2;re+=c*z[2*j]-s*z[2*j+1];im+=s*z[2*j]+c*z[2*j+1];}
 return re*re+im*im;
};
assert.ok(Math.abs(refined.peak-combined(offGrid)/(independent.reduce((s,z)=>s+combined(z),0)/independent.length))<1e-10);
console.log('PASS: Nelder–Mead improves the off-grid peak, recovers wrapped phases, and preserves the quiet-referenced statistic');
