"""Check coherent multi-chirp FFT accumulation against direct Eq.25 sums."""
from pathlib import Path
import h5py
import numpy as np
c=299792458.; f0=77e9; S=99.987e12; fs=12.5e6
x,v,y=.1,341.51,.15
u=2e-6+np.arange(225)/fs; uc=(u[0]+u[-1])/2; h=u-uc
nfft=32768
with h5py.File(Path(__file__).with_suffix('.h5'),'w') as out:
 out.attrs['model']='receive-time range, unit amplitude, coherent chirps, noiseless'; out.attrs['script']=Path(__file__).name
 for N in [1,2,4,8,16]:
  # Deliberately not an integer number of ADC samples; no compressed-clock FFT.
  starts=np.arange(N)*25.37e-6
  summed=0j; direct=0j; dat=[]; phases=[]; frequencies=[]
  for a in starts:
   t=a+u; R=np.hypot(x-v*t,y); tau=2*R/c
   phi=2*np.pi*(-f0*tau-S*u*tau+.5*S*tau*tau)
   tc=a+uc; Rc=np.hypot(x-v*tc,y); rd=-v*(x-v*tc)/Rc
   phic=2*np.pi*(-f0*2*Rc/c-S*uc*2*Rc/c+.5*S*(2*Rc/c)**2)
   omega=2*np.pi*(-2/c*((f0+S*uc)*rd+S*Rc)+4*S*Rc*rd/c**2)
   corr=np.exp(1j*(phi-phic-omega*h)); d=np.exp(1j*phi)
   k=int(np.round(omega/(2*np.pi)*nfft/fs)); fk=k*fs/nfft
   fft=np.fft.fft(d*np.conj(corr),nfft)[k%nfft]
   Y=fft*np.exp(-2j*np.pi*fk*(u[0]-uc))
   summed+=np.exp(-1j*phic)*Y
   exact=np.sum(d*np.conj(corr)*np.exp(-1j*omega*h))*np.exp(-1j*phic)
   direct+=np.sum(d*np.exp(-1j*phi))
   np.testing.assert_allclose(exact,np.sum(d*np.exp(-1j*phi)),atol=1e-8)
   dat.append(d);phases.append(phi);frequencies.append(omega/(2*np.pi))
  retention=abs(summed/direct)**2
  assert retention>.99
  g=out.create_group(str(N));g.attrs['power_retention']=retention
  for key,val in [('chirp_start_s',starts),('fast_time_s',u),('iq',dat),('phase',phases),('beat_frequency_hz',frequencies)]: g.create_dataset(key,data=val)
  print(N,'chirps: direct/FFT power retention',retention)
print('Multi-chirp exact-factorization and FFT checks passed')
