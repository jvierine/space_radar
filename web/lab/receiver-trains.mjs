// Mean-subtracted events; RAW common quiet trains for noise normalization. Overlap is allowed; references are correlated.
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
export function prepareReceiverTrains(receivers,{samples,rows,perFrame,start,pulses,bgStart,bgStop,noiseStart,noiseStop,fullBandwidth=false}) {
  if(receivers.length!==4 || receivers.some(z=>z.length!==2*samples*rows))throw Error('Need four synchronized receivers.');
  if(![start,pulses,bgStart,bgStop,noiseStart,noiseStop].every(Number.isInteger) || bgStart<0 || bgStart>=bgStop || bgStop>rows || noiseStart<0 || noiseStart>=noiseStop || noiseStop>rows)throw Error('Invalid event/background bounds.');
  const valid=commonValid(receivers,samples,rows);
  const trainValid=k=>k>=0&&k+pulses<=rows&&Math.floor(k/perFrame)===Math.floor((k+pulses-1)/perFrame)&&valid.slice(k,k+pulses).every(Boolean);
  if(!trainValid(start))throw Error('Selected train is not intact on all four receivers.');
  const meanCounts=[];
  const means=receivers.map(data=>{
    const ownValid=commonValid([data],samples,rows);
    const quiet=Array.from({length:bgStop-bgStart},(_,i)=>bgStart+i).filter(k=>ownValid[k]);
    if(!quiet.length)throw Error('Receiver has no intact background-mean chirps.');
    meanCounts.push(quiet.length);
    const mean=new Float64Array(2*samples);
    for(const k of quiet)for(let j=0;j<2*samples;j++)mean[j]+=data[2*k*samples+j]/quiet.length;
    return mean;
  });
  if(fullBandwidth){
    const quiet=Array.from({length:noiseStop-noiseStart},(_,i)=>noiseStart+i).filter(k=>valid[k]&&(k<start||k>=start+pulses));
    if(!quiet.length)throw Error('Blue background needs one intact chirp outside the analysis train.');
    // Remove the repeated wall waveform at each fast-time sample, separately per RX.
    // M/(M-1) corrects the variance lost by estimating that waveform from these samples.
    const ownCounts=[],powers=[];
    for(let a=0;a<4;a++){
      const own=commonValid([receivers[a]],samples,rows);
      const starts=Array.from({length:noiseStop-noiseStart},(_,i)=>noiseStart+i).filter(k=>own[k]&&(k<start||k>=start+pulses));
      if(!starts.length)throw Error(`RX${a} background has no intact chirp.`);
      ownCounts.push(starts.length);
      const mean=new Float64Array(2*samples);
      for(const k of starts)for(let j=0;j<2*samples;j++)mean[j]+=receivers[a][2*k*samples+j]/starts.length;
      let power=0;
      for(const k of starts)for(let j=0;j<2*samples;j++){
        const z=receivers[a][2*k*samples+j]-(starts.length>1?mean[j]:0);power+=z*z;
      }
      powers.push(power/(samples*Math.max(1,starts.length-1)));
    }
    const covariance=new Float64Array(32),count=quiet.length*samples;
    const commonMeans=receivers.map(z=>{
      const mean=new Float64Array(2*samples);
      if(quiet.length>1)for(const k of quiet)for(let j=0;j<2*samples;j++)mean[j]+=z[2*k*samples+j]/quiet.length;
      return mean;
    });
    const divisor=samples*Math.max(1,quiet.length-1);
    for(const k of quiet)for(let j=0;j<samples;j++)for(let a=0;a<4;a++)for(let b=0;b<4;b++){
      const t=2*(k*samples+j),ar=receivers[a][t]-commonMeans[a][2*j],ai=receivers[a][t+1]-commonMeans[a][2*j+1],br=receivers[b][t]-commonMeans[b][2*j],bi=receivers[b][t+1]-commonMeans[b][2*j+1],i=2*(a*4+b);
      covariance[i]+=(ar*br+ai*bi)/divisor;covariance[i+1]+=(ai*br-ar*bi)/divisor;
    }
    // Preserve a PSD correlation matrix while each diagonal uses its RX's own valid samples.
    const diagonal=powers.map((_,a)=>covariance[2*(4*a+a)]);
    for(let a=0;a<4;a++)for(let b=0;b<4;b++){
      const scale=diagonal[a]>0&&diagonal[b]>0?Math.sqrt(powers[a]*powers[b]/(diagonal[a]*diagonal[b])):0,i=2*(4*a+b);
      covariance[i]*=scale;covariance[i+1]*=scale;
      if(a===b)covariance[i]=powers[a];
    }
    const noiseCalibrated=ownCounts.every(n=>n>1),noiseMethod=noiseCalibrated?'background residual variance, unbiased mean correction':'single-chirp raw background upper bound';
    const width=2*samples*pulses,data=new Float32Array(width*8);
    for(let a=0;a<4;a++)for(let k=0;k<pulses;k++)for(let j=0;j<2*samples;j++)data[a*width+k*2*samples+j]=receivers[a][2*(start+k)*samples+j]-means[a][j];
    return {data,means,meanCounts,referenceStarts:quiet,quietCount:1,noisePower:Float32Array.from({length:4},(_,a)=>covariance[2*(a*4+a)]),noiseCovariance:covariance,noiseSamples:count,noiseSamplesPerRx:ownCounts.map(n=>n*samples),noiseCalibrated,noiseMethod};
  }
  const candidates=[];
  for(let k=noiseStart;k+pulses<=noiseStop;) {
    if(trainValid(k)&&(k+pulses<=start||k>=start+pulses)){candidates.push(k);k++;}else k++;
  }
  if(candidates.length<8)throw Error('Background must contain at least eight intact quiet windows outside the event (overlap allowed). Enlarge the blue window.');
  const count=Math.min(24,candidates.length);
  const referenceStarts=Array.from({length:count},(_,i)=>candidates[Math.floor(i*candidates.length/count)]);
  const width=2*samples*pulses,data=new Float32Array(width*4*(count+1));
  const put=(receiver,k,index,subtractMean)=>{
    for(let n=0;n<pulses;n++)for(let j=0;j<2*samples;j++)data[index*width+n*2*samples+j]=receivers[receiver][2*(k+n)*samples+j]-(subtractMean?means[receiver][j]:0);
  };
  for(let a=0;a<4;a++) {
    put(a,start,a,true);
    referenceStarts.forEach((k,i)=>put(a,k,4+a*count+i,false));
  }
  return {data,means,meanCounts,referenceStarts,quietCount:count};
}
