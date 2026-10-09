// Phase-only coherent receive beamforming, conditional on a trajectory template.
// Positive phase rotates the stored receiver voltage by exp(+i phase)/2.
export function phaseSearch(event, noise, steps = 10) {
  if (!Number.isInteger(steps) || steps < 2 || steps > 32 || event.length !== 8 || noise.length < 8)
    throw Error('Need four receivers, at least eight quiet trains, and 2–32 phase steps.');
  const scores = new Float64Array(steps ** 3), projection = new Float32Array(steps ** 2);
  projection.fill(-Infinity);
  const sum = (z, phases) => {
    let re=0,im=0;
    for(let a=0;a<4;a++) {
      const c=Math.cos(phases[a])/2,s=Math.sin(phases[a])/2;
      re+=c*z[2*a]-s*z[2*a+1]; im+=s*z[2*a]+c*z[2*a+1];
    }
    return re*re+im*im;
  };
  let peak=-Infinity, best=null;
  for(let i=0;i<steps;i++) for(let j=0;j<steps;j++) for(let k=0;k<steps;k++) {
    const phases=[0,i*2*Math.PI/steps,j*2*Math.PI/steps,k*2*Math.PI/steps];
    const reference=noise.reduce((total,z)=>total+sum(z,phases),0)/noise.length;
    const score=reference>0?sum(event,phases)/reference:NaN;
    const index=(i*steps+j)*steps+k; scores[index]=score;
    const db=10*Math.log10(score);
    if(Number.isFinite(db))projection[j*steps+i]=Math.max(projection[j*steps+i],db);
    if(Number.isFinite(score)&&score>peak){peak=score;best={indices:[i,j,k],phases,reference};}
  }
  if(!best)throw Error('No finite beam score: quiet matched power is zero.');
  const single=Array.from({length:4},(_,a)=>{
    const power=z=>z[2*a]**2+z[2*a+1]**2;
    return power(event)/(noise.reduce((s,z)=>s+power(z),0)/noise.length);
  });
  return {steps,peak, ...best, single, projection};
}

export function commonValid(receivers, samples, rows) {
  return Array.from({length:rows},(_,k)=>receivers.every(data=>{
    let run=0;
    for(let j=0;j<samples;j++) {
      const r=data[2*(k*samples+j)],i=data[2*(k*samples+j)+1];
      if(!Number.isFinite(r)||!Number.isFinite(i))return false;
      run=(r===0&&i===0)?run+1:0;
      if(run>=8)return false;
    }
    return true;
  }));
}

export function projectReceivers(receivers, q, {samples, rows, perFrame, start, pulses, bgStart, bgStop, noiseStart, noiseStop, conjugated}) {
  const valid=commonValid(receivers,samples,rows);
  const trainValid=k=>k>=0&&k+pulses<=rows&&Math.floor(k/perFrame)===Math.floor((k+pulses-1)/perFrame)&&valid.slice(k,k+pulses).every(Boolean);
  if(!trainValid(start))throw Error('Selected train is not intact on all four receivers.');
  const quiet=Array.from({length:bgStop-bgStart},(_,i)=>bgStart+i).filter(k=>valid[k]);
  if(!quiet.length)throw Error('No common intact quiet-mean chirps across four receivers.');
  const means=receivers.map(data=>{
    const mean=new Float64Array(2*samples);
    for(const k of quiet)for(let j=0;j<2*samples;j++)mean[j]+=data[2*k*samples+j]/quiet.length;
    return mean;
  });
  const project=k=>{
    const out=new Float64Array(8);
    for(let a=0;a<4;a++)for(let n=0;n<pulses;n++)for(let j=0;j<samples;j++){
      const ti=2*(n*samples+j),zi=2*((k+n)*samples+j);
      const qr=q[ti],qi=(conjugated?-1:1)*q[ti+1];
      const zr=receivers[a][zi]-means[a][2*j],zii=receivers[a][zi+1]-means[a][2*j+1];
      out[2*a]+=qr*zr+qi*zii;out[2*a+1]+=qr*zii-qi*zr;
    }
    return out;
  };
  const candidates=[];
  for(let k=noiseStart;k+pulses<=noiseStop;){
    // Independent reference blocks, outside event AND background-training data.
    if(trainValid(k)&&(k+pulses<=start||k>=start+pulses)&&(k+pulses<=bgStart||k>=bgStop)){candidates.push(k);k+=pulses;}else k++;
  }
  if(candidates.length<8)throw Error('Need at least eight common intact noise trains outside the event and mean interval.');
  const count=Math.min(24,candidates.length);
  const referenceStarts=Array.from({length:count},(_,i)=>candidates[Math.floor(i*candidates.length/count)]);
  return {event:project(start),noise:referenceStarts.map(project),referenceStarts,meanCount:quiet.length};
}
