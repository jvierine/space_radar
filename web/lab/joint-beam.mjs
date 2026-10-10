import {noiseMetric,metricProduct,metricEnergy,fittedNoiseBias} from './noise-metric.mjs';
// Variable-projection least squares: A=s/Q. Common temporal template,
// constant receive covariance, temporally white noise on acquired samples.
export function jointBeam({event,noiseCovariance,templateEnergy=1}){
 const {precision,regularization}=noiseMetric(noiseCovariance);
 const glsScore=metricEnergy(precision,event),bias=fittedNoiseBias(precision,noiseCovariance);
 const signalSnr=Math.max(0,glsScore-bias);
 const unnormalized=metricProduct(precision,event),norm=Math.hypot(...unnormalized);
 if(!(norm>0))throw Error('Zero fitted echo amplitude; no beam direction can be estimated.');
 const weights=Float64Array.from(unnormalized,v=>v/norm);
 const reference=metricEnergy(noiseCovariance,weights);
 const amplitudes=Float64Array.from(event,v=>v/templateEnergy);
 const phases=Array.from({length:4},(_,b)=>Math.atan2(event[2*b+1],event[2*b]));
 const echoPhases=phases.map(p=>Math.atan2(Math.sin(p-phases[0]),Math.cos(p-phases[0])));
 const channelPower=Array.from({length:4},(_,b)=>({observed:event[2*b]**2+event[2*b+1]**2,quiet:noiseCovariance[2*(4*b+b)]}));
 const single=channelPower.map(p=>p.observed/p.quiet);
 // Equal calibrated antenna response measured in units of each RX noise RMS.
 const response=new Float64Array(8);
 for(let b=0;b<4;b++){const rms=Math.sqrt(noiseCovariance[2*(4*b+b)]);response[2*b]=rms*Math.cos(phases[b]);response[2*b+1]=rms*Math.sin(phases[b]);}
 // Actual implemented beam response/noise, including any covariance ridge.
 let re=0,im=0;for(let b=0;b<4;b++){re+=weights[2*b]*response[2*b]+weights[2*b+1]*response[2*b+1];im+=weights[2*b]*response[2*b+1]-weights[2*b+1]*response[2*b];}
 const rcsGain=(re*re+im*im)/reference;
 const growth=Array.from({length:4},(_,i)=>{const count=i+1,m=noiseMetric(noiseCovariance,count),s=event.slice(0,2*count),score=metricEnergy(m.precision,s),noise=fittedNoiseBias(m.precision,noiseCovariance,count);return {count,score:1+Math.max(score-noise,0),glsScore:score,fittedNoiseBias:noise};});
 return {peak:1+signalSnr,glsScore,fittedNoiseBias:bias,signalSnr,regularization,weights,amplitudes,echoPhases,phases:echoPhases,reference,single,channelPower,growth,rcsGain,noiseCovariance,noiseModel:'joint GLS with full background receive covariance; white temporal noise',objective:'min sum (z-q A)^H C^-1 (z-q A); A solved analytically, trajectory refined by bounded Nelder-Mead'};
}
