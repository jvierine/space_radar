import assert from 'node:assert/strict';
import {prepareReceiverTrains} from '../web/lab/receiver-trains.mjs';
const samples=64,rows=80,M=29,options={samples,rows,perFrame:80,start:40,pulses:8,bgStart:0,bgStop:M,noiseStart:0,noiseStop:M,fullBandwidth:true};
let seed=72;const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296-.5;};
const noise=Array.from({length:4},(_,a)=>Float32Array.from({length:2*samples*rows},()=>random()*(a+1)*8));
const clutter=noise.map((data,a)=>Float32Array.from(data,(value,i)=>value+(i%2?Math.sin(.2*Math.floor(i/2)%samples):Math.cos(.1*Math.floor(i/2)%samples))*0));
// Exactly repeat a strong arbitrary complex wall waveform at each fast-time sample.
for(let a=0;a<4;a++)for(let k=0;k<rows;k++)for(let j=0;j<samples;j++){
 clutter[a][2*(k*samples+j)]+=10000*(a+1)*Math.cos(.31*j);
 clutter[a][2*(k*samples+j)+1]+=10000*(a+1)*Math.sin(.31*j);
}
const pure=prepareReceiverTrains(noise,options),wall=prepareReceiverTrains(clutter,options);
for(let a=0;a<4;a++){
 let sum=0;const mean=new Float64Array(2*samples);
 for(let k=0;k<M;k++)for(let j=0;j<2*samples;j++)mean[j]+=clutter[a][2*k*samples+j]/M;
 for(let k=0;k<M;k++)for(let j=0;j<2*samples;j++)sum+=(clutter[a][2*k*samples+j]-mean[j])**2;
 assert(Math.abs(wall.noisePower[a]/(sum/(samples*(M-1)))-1)<1e-6);
 assert(Math.abs(wall.noisePower[a]/pure.noisePower[a]-1)<.001,'Stationary walls do not become thermal noise');
}
assert.equal(wall.noiseCalibrated,true);assert.deepEqual(wall.noiseSamplesPerRx,[M*samples,M*samples,M*samples,M*samples]);
const changed=noise.map(z=>z.slice());for(let j=0;j<2*samples*M;j++)changed[2][j]*=5;
const own=prepareReceiverTrains(changed,options);for(let a=0;a<4;a++)assert(Math.abs(own.noisePower[a]/pure.noisePower[a]-(a===2?25:1))<1e-5);
const single=prepareReceiverTrains(clutter,{...options,bgStop:1,noiseStop:1});assert.equal(single.noiseCalibrated,false);assert(single.noisePower.every(v=>v>0));
console.log('PASS: strong repeated walls removed, unbiased residual variance, all selected samples, independent RX noise, one-chirp upper-bound flag');
