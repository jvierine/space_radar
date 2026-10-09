"""Count the shipped radial bank and sparse FFT butterflies; no timing claims."""
from pathlib import Path
import json, math
import h5py
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
p=json.loads((ROOT/'web/lab/datasets/test63/metadata.json').read_text())['parameters']
c=299792458.; n=225; period=25.37e-6; fs=p['fs']; slope=p['freq_slope']; f=p['f_start']+slope*(p['T_adc']+(n-1)/(2*fs))
epsilon=math.acos(math.sqrt(.95)); nf=2048; rlo,rhi=.001,3.; vlo,vhi=0.,600.; alo,ahi=0.,1e6
rows=[]
for nc in [1,2,4,8,16,32]:
 ns=1 if nc==1 else 2**math.ceil(math.log2(8*nc))
 h=(nc-1)*period/2+(n-1)/(2*fs);eta=(n-1)/(2*fs);k=4*math.pi/c;l=4*math.pi*slope/c**2
 da=.5*(k*abs(f)+2*abs(l)*rhi)*h*h+.5*k*abs(slope)*eta*h*h+abs(l)*(vhi*h**3+.5*ahi*h**4)
 dv=k*abs(slope)*eta*h+abs(l)*(2*vhi*h*h+ahi*h**3)
 dr=abs(l)*ahi*h*h*(rhi-rlo)/2
 astep=2*epsilon/max(da,1e-30);vg=2*(epsilon-dr)/max(dv,1e-30)
 rstep=c*fs/(2*abs(slope)*nf);vstep=c*fs/(2*abs(f)*nf) if ns==1 else c/(2*abs(f)*ns*period)
 nr=math.ceil((rhi-rlo)/rstep)+1;nv=math.ceil((vhi-vlo)/vstep)+1;na=math.ceil((ahi-alo)/astep)+1
 rs=np.linspace(rlo,rhi,nr);vs=np.linspace(vlo,vhi,nv)
 groups=[];first=0
 while first<nv:
  last=first+1
  while last<nv and vs[last]-vs[first]<=vg:last+=1
  ff=(-2*(p['f_start']+slope*(p['T_adc']+eta))/c+4*slope*rs[:,None]/c**2)*vs[None,first:last]-2*slope*rs[:,None]/c
  # Count before invalid-positive-range masking: a conservative work bound.
  bins=np.remainder(np.rint(ff/fs*nf).astype(np.int64),nf)
  groups.append(len(np.unique(bins)));first=last
 G=na*len(groups);nodes=nr*nv*na
 direct=8*nodes*n*nc
 butterflies=8*na*sum(nc*nf/2*math.log2(nf)+cols*ns/2*math.log2(ns) for cols in groups)
 corrections=8*G*n*nc;reads=8*nodes
 rows.append([nc,nr,nv,na,G,nodes,ns,direct,butterflies,corrections,reads,nodes*4,nc*n/fs,dr<epsilon and nodes<=32000000])
 print(nc,nr,nv,na,G,nodes, f'direct {direct:.3g}, FFT butterflies {butterflies:.3g}', 'accepted' if rows[-1][-1] else 'over cap')
output=ROOT/'manim/assets/radial_complexity.h5'
with h5py.File(output,'w') as f:
 f['counts']=np.array(rows,dtype=np.float64)
 f.attrs['columns']='chirps,nr,nv,na,groups,nodes,slow_fft,direct_complex_MAC,fft_butterflies,correction_complex_mult,score_reads,cube_bytes,T_coh_s,accepted'
 f.attrs['source']='tools/radial_complexity.py; lab-core/src/radial_search.rs'
 f.attrs['bounds']='r0 .001..3 m; v0 0..600 m/s; a0 0..1e6 m/s2; loss .05'
np.set_printoptions(suppress=True)
