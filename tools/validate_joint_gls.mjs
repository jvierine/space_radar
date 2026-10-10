import assert from 'node:assert/strict';
import {noiseMetric,metricProduct,metricEnergy} from '../web/lab/noise-metric.mjs';
import {jointBeam} from '../web/lab/joint-beam.mjs';
const C=new Float64Array(32);const diag=[4,9,16,25];
for(let i=0;i<4;i++)C[2*(4*i+i)]=diag[i];C[2]=1;C[3]=1;C[8]=1;C[9]=-1;
const {precision}=noiseMetric(C);
const expected=new Float64Array(32);expected[0]=9/34;expected[2]=-1/34;expected[3]=-1/34;expected[8]=-1/34;expected[9]=1/34;expected[10]=4/34;expected[20]=1/16;expected[30]=1/25;
precision.forEach((v,i)=>assert.ok(Math.abs(v-expected[i])<1e-13));
const s=Float64Array.of(12,3,-2,8,5,-7,8,1),Q=7;
const beam=jointBeam({event:s,noiseCovariance:Float64Array.from(C,v=>v*Q),templateEnergy:Q});
assert.ok(Math.abs(beam.fittedNoiseBias-4)<1e-12);
assert.ok(Math.abs(beam.glsScore-metricEnergy(precision,s)/Q)<1e-12);
let re=0,im=0;for(let i=0;i<4;i++){re+=beam.weights[2*i]*s[2*i]+beam.weights[2*i+1]*s[2*i+1];im+=beam.weights[2*i]*s[2*i+1]-beam.weights[2*i+1]*s[2*i];}
assert.ok(Math.abs((re*re+im*im)/beam.reference-beam.glsScore)<1e-11);
// Independent direct residual evaluation: J(A)=constant-score at A=s/Q.
const q=[1,0,0,1,-1,0],zs=[[3,1,1,-2,4,2,-1,3],[0,2,3,1,2,-3,1,1],[-2,1,0,1,1,2,2,1]];
const sums=new Float64Array(8);for(let j=0;j<3;j++)for(let b=0;b<4;b++){sums[2*b]+=q[2*j]*zs[j][2*b]+q[2*j+1]*zs[j][2*b+1];sums[2*b+1]+=q[2*j]*zs[j][2*b+1]-q[2*j+1]*zs[j][2*b];}
const A=Float64Array.from(sums,v=>v/3),constant=zs.reduce((v,z)=>v+metricEnergy(precision,z),0);
const residual=zs.reduce((v,z,j)=>v+metricEnergy(precision,z.map((x,i)=>{const b=i-i%2;return x-(i%2===0?q[2*j]*A[b]-q[2*j+1]*A[b+1]:q[2*j]*A[b+1]+q[2*j+1]*A[b]);})),0);
assert.ok(Math.abs(residual-(constant-metricEnergy(precision,sums)/3))<1e-12);
const equalC=new Float64Array(32);for(let b=0;b<4;b++)equalC[2*(4*b+b)]=1;
const equal=jointBeam({event:Float64Array.of(10,0,0,10,-10,0,0,-10),noiseCovariance:equalC});assert.ok(Math.abs(equal.rcsGain-4)<1e-12);
assert.throws(()=>noiseMetric(Float64Array.from(C,(v,i)=>i===3?NaN:v)));
console.log('PASS: analytic covariance inverse, direct weighted residual identity, complex beam/noise propagation, fitted-amplitude bias and equal-RX calibration');
