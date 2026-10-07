// Equation 25. Signed convention: received RF times conjugate of transmit LO.
export const C = 299792458;
export const PROFILES = {
  '6607': {name:'6607 · 79.5 GHz · 60 MHz/µs', f0:79.5e9, slope:60e12, fs:12.5e6, samples:128, adcStart:0, tail:0.37e-6, idle:3.099999904632568e-6, speed:6285},
  '6606': {name:'6606 · 78 GHz · 100 MHz/µs', f0:78e9, slope:100e12, fs:12.5e6, samples:128, adcStart:0, tail:0.37e-6, idle:3.099999904632568e-6, speed:6700}
};
export function timing(p) {
  const live=p.samples/p.fs, ramp=p.adcStart+live+p.tail;
  return {live,ramp,period:ramp+p.idle};
}
export function state(t,u,p,g) {
  const x=g.v*t-g.x0, rr=x*x+g.y0*g.y0;
  // Stable positive solution for one-way retarded flight time.
  const d=rr/(Math.sqrt((x*g.v)**2+(C*C-g.v*g.v)*rr)+x*g.v);
  const tau=2*d, xs=x-g.v*d, range=Math.hypot(xs,g.y0), vr=g.v*xs/range;
  const tauDot=2*vr/(C+vr), doppler=-p.f0*tauDot, beat=-p.slope*tau;
  const motion=-p.slope*(u-tau)*tauDot, frequency=beat+doppler+motion;
  const phase=2*Math.PI*(-p.f0*tau-p.slope*u*tau+.5*p.slope*tau*tau);
  return {tau,range,vr,doppler,beat,motion,frequency,phase,amplitude:1/(range*range),sameRamp:u>=tau};
}
export function alias(f,fs) {return ((f+fs/2)%fs+fs)%fs-fs/2;}
export function voltage(s,mode='ideal') {
  if(!s.sameRamp) return [0,0];
  let a=s.amplitude, re=Math.cos(s.phase), im=Math.sin(s.phase);
  if(mode==='receiver') {
    // Teaching approximation: quasistatic two first-order analogue HPFs,
    // followed by ideal range-side Complex-1x visibility mask.
    if(s.frequency>=0 || s.frequency<=-10e6) return [0,0];
    for(const fc of [175e3,350e3]) {
      const den=fc*fc+s.frequency*s.frequency, hr=s.frequency*s.frequency/den, hi=fc*s.frequency/den;
      const r=re*hr-im*hi; im=re*hi+im*hr;re=r;
    }
  }
  return [a*re,a*im];
}
export function simulate(p,g,n,mode='ideal') {
  const clock=timing(p), records=[], dense=[], points=[];
  for(let k=0;k<n;k++) {
    const record=[];
    for(let j=0;j<p.samples;j++) {
      const u=p.adcStart+j/p.fs, t=k*clock.period+u, s=state(t,u,p,g), [i,q]=voltage(s,mode);
      const point={t,u,k,j,...s,i,q,aliased:alias(s.frequency,p.fs)};
      record.push(point);points.push(point);
    }
    records.push(record);
    // Dense analytic curves preserve high Doppler even above ADC Nyquist.
    const analytic=[];
    for(let j=0;j<=256;j++) {
      const u=clock.ramp*j/256,t=k*clock.period+u;
      analytic.push({t,u,k,...state(t,u,p,g),tx:p.f0+p.slope*u});
    }
    dense.push(analytic);
  }
  return {p,g,n,mode,clock,records,points,dense,duration:n*clock.period,span:(n-1)*clock.period+p.adcStart+clock.live};
}
export function fft(re,im) {
  const n=re.length;
  for(let i=1,j=0;i<n;i++) {let bit=n>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j){[re[i],re[j]]=[re[j],re[i]];[im[i],im[j]]=[im[j],im[i]];}}
  for(let len=2;len<=n;len<<=1) {
    const angle=-2*Math.PI/len,wr0=Math.cos(angle),wi0=Math.sin(angle);
    for(let base=0;base<n;base+=len) {
      let wr=1,wi=0;
      for(let j=0;j<len/2;j++) {
        const a=base+j,b=a+len/2, tr=wr*re[b]-wi*im[b],ti=wr*im[b]+wi*re[b];
        re[b]=re[a]-tr;im[b]=im[a]-ti;re[a]+=tr;im[a]+=ti;
        const w=wr*wr0-wi*wi0;wi=wr*wi0+wi*wr0;wr=w;
      }
    }
  }
}
export function spectrogram(sim,size=64,hop=16,nfft=256) {
  const frames=[];let peak=0;
  for(const record of sim.records) {
    for(let start=0;start+size<=record.length;start+=hop) {
      const re=new Float64Array(nfft),im=new Float64Array(nfft);
      for(let j=0;j<size;j++) {const w=.5-.5*Math.cos(2*Math.PI*j/(size-1));re[j]=record[start+j].i*w;im[j]=record[start+j].q*w;}
      fft(re,im);const power=new Float64Array(nfft);
      for(let b=0;b<nfft;b++) {const ix=(b+nfft/2)%nfft;power[b]=re[ix]*re[ix]+im[ix]*im[ix];peak=Math.max(peak,power[b]);}
      frames.push({t:record[start+(size>>1)].t,k:record[0].k,power});
    }
  }
  return {frames,size,hop,nfft,peak,fs:sim.p.fs,resolution:sim.p.fs/size};
}
