import {metricEnergy} from './noise-metric.mjs';
import {snrDb} from './rcs.mjs';

// Freeze the fitted spatial beam and motion; change only midpoint range.
export function coherentRangeProfile({receivers,means,start,pulses,samples,weights,rawCovariance,ranges,templateAt,conjugated=false}) {
  const voltage=new Float64Array(2*pulses*samples);
  for(let n=0;n<pulses;n++)for(let j=0;j<samples;j++)for(let b=0;b<4;b++){
    const src=2*((start+n)*samples+j),dst=2*(n*samples+j);
    const re=receivers[b][src]-means[b][2*j],im=receivers[b][src+1]-means[b][2*j+1];
    voltage[dst]+=weights[2*b]*re+weights[2*b+1]*im;
    voltage[dst+1]+=weights[2*b]*im-weights[2*b+1]*re;
  }
  const noisePerSample=metricEnergy(rawCovariance,weights);
  if(!(noisePerSample>0))throw Error('Coherent range map needs positive beam noise power.');
  const observedPower=new Float64Array(ranges.length),noisePower=new Float64Array(ranges.length),score=new Float64Array(ranges.length),snr=new Float32Array(ranges.length);
  for(let r=0;r<ranges.length;r++){
    const q=templateAt(ranges[r]);
    if(q.length!==voltage.length)throw Error('Range template length mismatch.');
    let re=0,im=0,energy=0;
    for(let j=0;j<q.length;j+=2){
      const qr=q[j],qi=(conjugated?-1:1)*q[j+1];
      re+=qr*voltage[j]+qi*voltage[j+1];im+=qr*voltage[j+1]-qi*voltage[j];energy+=qr*qr+qi*qi;
    }
    if(!(energy>0)){observedPower[r]=noisePower[r]=score[r]=snr[r]=NaN;continue;}
    observedPower[r]=(re*re+im*im)/(energy*energy);
    noisePower[r]=noisePerSample/energy;
    score[r]=observedPower[r]/noisePower[r];snr[r]=snrDb(score[r]);
  }
  return {ranges,observedPower,noisePower,score,snrDb:snr,weights:Float64Array.from(weights),normalization:'|q^H (w^H z)|^2 / [(q^H q)(w^H C w)]; displayed SNR subtracts one expected noise unit and floors at 0 dB. Peak beam weights, velocity and acceleration held fixed across range; selection bias remains.'};
}
