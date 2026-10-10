// Complex Hermitian receive-noise metric. Temporal noise is assumed white.
const mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]];
const sub=(a,b)=>[a[0]-b[0],a[1]-b[1]];
const conj=a=>[a[0],-a[1]];
export function noiseMetric(covariance,count=4){
 if(covariance.length<32||Array.from(covariance).some(x=>!Number.isFinite(x)))throw Error('Finite four-RX covariance required.');
 const d=Array.from({length:count},(_,i)=>Math.sqrt(covariance[2*(4*i+i)]));
 if(d.some(x=>!(x>0&&Number.isFinite(x))))throw Error('Positive per-RX background noise is required.');
 const R=Array.from({length:count},(_,i)=>Array.from({length:count},(_,j)=>{
  const a=2*(4*i+j),b=2*(4*j+i);
  return [(covariance[a]+covariance[b])/(2*d[i]*d[j]),(covariance[a+1]-covariance[b+1])/(2*d[i]*d[j])];
 }));
 for(const ridge of [0,1e-8,1e-6,1e-4,.01]){
  const L=Array.from({length:count},()=>Array.from({length:count},()=>[0,0]));let valid=true;
  for(let i=0;i<count&&valid;i++)for(let j=0;j<=i;j++){
   let z=[R[i][j][0]+(i===j?ridge:0),R[i][j][1]];
   for(let k=0;k<j;k++)z=sub(z,mul(L[i][k],conj(L[j][k])));
   if(i===j){if(!(z[0]>1e-9)){valid=false;break;}L[i][j]=[Math.sqrt(z[0]),0];}
   else L[i][j]=z.map(v=>v/L[j][j][0]);
  }
  if(!valid)continue;
  const precision=new Float64Array(2*count*count);
  for(let column=0;column<count;column++){
   const y=[],x=Array.from({length:count},()=>[0,0]);
   for(let i=0;i<count;i++){let z=[i===column?1:0,0];for(let j=0;j<i;j++)z=sub(z,mul(L[i][j],y[j]));y[i]=z.map(v=>v/L[i][i][0]);}
   for(let i=count-1;i>=0;i--){let z=y[i];for(let j=i+1;j<count;j++)z=sub(z,mul(conj(L[j][i]),x[j]));x[i]=z.map(v=>v/L[i][i][0]);}
   for(let i=0;i<count;i++){const offset=2*(i*count+column);precision[offset]=x[i][0]/(d[i]*d[column]);precision[offset+1]=x[i][1]/(d[i]*d[column]);}
  }
  return {precision,regularization:ridge,count};
 }
 throw Error('Background covariance is not positive semidefinite.');
}
export function metricProduct(precision,vector,count=vector.length/2){
 const out=new Float64Array(2*count);
 for(let i=0;i<count;i++)for(let j=0;j<count;j++){
  const k=2*(i*count+j),p=precision[k],q=precision[k+1],r=vector[2*j],s=vector[2*j+1];
  out[2*i]+=p*r-q*s;out[2*i+1]+=p*s+q*r;
 }
 return out;
}
export function metricEnergy(precision,vector){
 const out=metricProduct(precision,vector);
 return vector.reduce((sum,v,i)=>sum+v*out[i],0);
}
export function fittedNoiseBias(precision,covariance,count=4){
 let trace=0;for(let i=0;i<count;i++)for(let j=0;j<count;j++){
  const a=2*(i*count+j),b=2*(4*j+i);trace+=precision[a]*covariance[b]-precision[a+1]*covariance[b+1];
 }
 return trace;
}
