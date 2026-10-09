const C=299792458,KB=1.380649e-23;

// PEC-sphere monostatic Mie series, using Riccati-Bessel coefficients.
// Same convention as capture_analysis/analyze_captures.py:sphere.
export function sphereRcs(frequency,diameter) {
  if(!(frequency>0 && diameter>0))return NaN;
  const k=2*Math.PI*frequency/C,x=k*diameter/2;
  const nmax=Math.ceil(x+4*Math.cbrt(x)+12),m=nmax+50;
  const j=new Float64Array(m+2);j[m]=1;
  for(let n=m;n>0;n--) {
    j[n-1]=(2*n+1)/x*j[n]-j[n+1];
    if(Math.abs(j[n-1])>1e100)for(let i=n-1;i<=m+1;i++)j[i]*=1e-100;
  }
  const j0=Math.sin(x)/x,j1=(j0-Math.cos(x))/x;
  const scale=Math.abs(j0)>Math.abs(j1)?j0/j[0]:j1/j[1];
  for(let n=0;n<=nmax;n++)j[n]*=scale;
  let previousY=-Math.cos(x)/x,y=previousY/x-Math.sin(x)/x,re=0,im=0;
  const ratio=(a,b)=>{const norm=a*a+b*b;return [-a*a/norm,a*b/norm];};
  for(let n=1;n<=nmax;n++) {
    const a=ratio(x*j[n-1]-n*j[n],x*previousY-n*y),b=ratio(j[n],y);
    const weight=(2*n+1)*(n%2?-1:1);
    re+=weight*(a[0]-b[0]);im+=weight*(a[1]-b[1]);
    const nextY=(2*n+1)/x*y-previousY;previousY=y;y=nextY;
  }
  return Math.PI/(k*k)*(re*re+im*im);
}

export function diameterRoots(sigma,frequency,minimum,maximum) {
  if(!(sigma>0 && minimum>0 && maximum>=minimum))return [];
  const roots=[],intervals=4096;
  const value=d=>sphereRcs(frequency,d)-sigma;
  let a=minimum,fa=value(a);
  if(fa===0)roots.push(a);
  for(let i=1;i<=intervals;i++) {
    const b=minimum*Math.pow(maximum/minimum,i/intervals),fb=value(b);
    if(fa*fb<0) {
      let lo=a,hi=b,flo=fa;
      for(let j=0;j<48;j++) {
        const middle=(lo+hi)/2,fm=value(middle);
        if(flo*fm<=0)hi=middle;else {lo=middle;flo=fm;}
      }
      roots.push((lo+hi)/2);
    } else if(fb===0)roots.push(b);
    a=b;fa=fb;
  }
  return roots.filter((d,i)=>i===0 || Math.abs(d-roots[i-1])>1e-10);
}

export function estimateRcs(ratio,{temperature,range,frequency,time,txPowerDbm,txGainDbi,rxGainDbi,lossDb},receivers=1) {
  const signalSnr=Math.max(ratio-1,0);
  const power=10**((txPowerDbm-30)/10),gain=10**((txGainDbi+rxGainDbi-lossDb)/10);
  return signalSnr*KB*temperature/time*(4*Math.PI)**3*range**4/(power*gain*(C/frequency)**2*receivers);
}
