import assert from 'node:assert/strict';
import fs from 'node:fs';
import vm from 'node:vm';
import {coherentRangeProfile} from '../web/lab/coherent-range.mjs';
const N=64,phases=[0,.4,-1,2],amplitude=3;
const templateAt=bin=>Float64Array.from({length:2*N},(_,i)=>i%2?Math.sin(2*Math.PI*bin*Math.floor(i/2)/N):Math.cos(2*Math.PI*bin*Math.floor(i/2)/N));
const covariance=new Float64Array(32);for(let b=0;b<4;b++)covariance[2*(4*b+b)]=2;
for(const conjugated of [false,true]){
 const receivers=phases.map(p=>Float64Array.from({length:2*N},(_,i)=>{
   const angle=(conjugated?-1:1)*2*Math.PI*8*Math.floor(i/2)/N+p;
   return amplitude*(i%2?Math.sin(angle):Math.cos(angle));
 }));
 const weights=Float64Array.from(phases.flatMap(p=>[Math.cos(p)/2,Math.sin(p)/2]));
 const p=coherentRangeProfile({receivers,means:Array.from({length:4},()=>new Float64Array(2*N)),start:0,pulses:1,samples:N,weights,rawCovariance:covariance,ranges:[8,12],templateAt,conjugated});
 assert(Math.abs(p.observedPower[0]-4*amplitude**2)<1e-10);
 assert(Math.abs(p.noisePower[0]-2/N)<1e-12);
 assert(Math.abs(p.score[0]-4*N*amplitude**2/2)<1e-8,'Four coherent antennas give four times the single-RX SNR');
 assert(p.observedPower[1]<1e-20,'Wrong range template cancels coherently');
 assert.equal(p.snrDb[1],0);
 assert.deepEqual([...p.weights],[...weights],'Peak phasing stays fixed for all ranges');
}
console.log('PASS: coherent four-RX gain, complex convention, range sidelobe rejection, propagated noise, fixed peak weights');
const source=fs.readFileSync('web/lab/app.mjs','utf8');
const context={period:.01,meta:{parameters:{T_adc:0}},clock:k=>(k*.01).toFixed(6),$:()=>({hidden:true}),coherentMap:{set(image,columns,rows,config){this.output={image,columns,rows,config};}}};
vm.createContext(context);vm.runInContext(source.slice(source.indexOf('function drawCoherentMap('),source.indexOf('function drawCoherentProfile(')),context);
context.points=[{start:10,midpoint:0,rangeProfile:{ranges:[.2,.3],snrDb:[1,2]}},{start:12,midpoint:0,rangeProfile:{ranges:[.2,.3],snrDb:[3,4]}}];
vm.runInContext('drawCoherentMap(points)',context);
const map=context.coherentMap.output;assert.equal(map.columns,3);assert.equal(map.rows,2);
assert.deepEqual([...map.image],[1,NaN,3,2,NaN,4],'Range is vertical, time horizontal, skipped trains remain blank');
assert.equal(map.config.colorLabel,'Beamformed SNR (dB)');
console.log('PASS: shipped GUI time/range orientation, blank missing windows and SNR color scale');
