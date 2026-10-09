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
    const covariance=new Float64Array(32),count=quiet.length*samples;
    for(const k of quiet)for(let j=0;j<samples;j++)for(let a=0;a<4;a++)for(let b=0;b<4;b++){
      const t=2*(k*samples+j),ar=receivers[a][t],ai=receivers[a][t+1],br=receivers[b][t],bi=receivers[b][t+1],i=2*(a*4+b);
      covariance[i]+=(ar*br+ai*bi)/count;covariance[i+1]+=(ai*br-ar*bi)/count;
    }
    const width=2*samples*pulses,data=new Float32Array(width*8);
    for(let a=0;a<4;a++)for(let k=0;k<pulses;k++)for(let j=0;j<2*samples;j++)data[a*width+k*2*samples+j]=receivers[a][2*(start+k)*samples+j]-means[a][j];
    return {data,means,meanCounts,referenceStarts:quiet,quietCount:1,noisePower:Float32Array.from({length:4},(_,a)=>covariance[2*(a*4+a)]),noiseCovariance:covariance,noiseSamples:count};
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
