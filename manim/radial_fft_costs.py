"""Analytical FFT-bank work model; no measured throughput. By Juha Vierinen."""
from pathlib import Path
import json, math
import numpy as np
import h5py

def costs(loss=1-10**(-2/10),save=True):
    p=json.loads((Path(__file__).parent.parent/'web/lab/datasets/test63/metadata.json').read_text())['parameters']
    fs=p['fs'];gamma=p['freq_slope'];c=299792458.;n=225;period=25.37e-6
    nf=2048;eta=(n-1)/(2*fs);fc=p['f_start']+gamma*(p['T_adc']+eta)
    k=4*math.pi/c;l=4*math.pi*gamma/c**2
    rows=[]
    for K in [1,2,4,8,16,32]:
        u=np.arange(n)/fs-eta
        h=(np.arange(K)[:,None]-(K-1)/2)*period+u[None,:]
        basis=[h,h*h,h**3,h**4,u*h,u*h*h]
        sd=[np.std(x) for x in basis];su=np.sqrt((n*n-1)/12)/fs
        da=.5*(k*abs(fc)+2*abs(l)*3)*sd[1]+.5*k*abs(gamma)*sd[5]+abs(l)*(600*sd[2]+.5e6*sd[3])
        dv=k*abs(gamma)*sd[4]+abs(l)*(1200*sd[1]+1e6*sd[2])
        dr=abs(l)*1e6*sd[1]*2.999/2
        gr=k*abs(gamma)*su+2*abs(l)*600*sd[0]+abs(l)*1e6*sd[1]
        gv=(k*abs(fc)+2*abs(l)*3)*sd[0]+dv
        nf=2048;ns=1 if K==1 else 8*K
        while True:
            rs=c*fs/(2*abs(gamma)*nf)
            vs=c*fs/(2*abs(fc)*nf) if ns==1 else c/(2*abs(fc)*ns*period)
            sf=math.pi*np.sqrt((n*n-1)/12)/nf
            ss=0 if ns==1 else math.pi*np.sqrt((K*K-1)/12)/ns
            lookup=math.hypot(sf,ss);grid=.5*(gr*rs+gv*vs)
            if lookup+grid<.6*math.sqrt(loss):
                epsilon=math.sqrt(loss)-lookup-grid;break
            nf*=2
            if ns>1:ns*=2
        astep=epsilon/max(da,1e-30);vgroup=(epsilon-2*dr)/max(dv,1e-30)
        assert vgroup>0
        nr=math.ceil(2.999/(c*fs/(2*gamma*nf)))+1
        nv=math.ceil(600/(c*fs/(2*fc*nf) if K==1 else c/(2*fc*ns*period)))+1
        na=math.ceil(1e6/astep)+1
        ranges=np.linspace(.001,3,nr);velocities=np.linspace(0,600,nv)
        first=0;Ms=[]
        while first<nv:
            last=first+1
            while last<nv and velocities[last]-velocities[first]<=vgroup:last+=1
            r=ranges[:,None];v=velocities[None,first:last]
            slow=(-2*fc/c+4*gamma*r/c**2)*v
            fast=slow-2*gamma*r/c
            # Union including invalid near-zero-range trajectories: upper column count.
            bins=np.rint(fast/fs*nf).astype(np.int64)%nf
            Ms.append(len(np.unique(bins)));first=last
        G=na*len(Ms);sumM=na*sum(Ms)
        fft_units=G*K*nf*math.log2(nf)+sumM*ns*math.log2(ns)
        corrections=4*G*K*n*6
        scoring=16*nr*nv*na
        real_ops=4*5*fft_units+corrections+scoring
        rows.append([K,na,G,nr*nv*na,real_ops,real_ops/period,sumM])
    path=Path(__file__).parent/'assets/radial_fft_costs.h5'
    if not save:return rows
    with h5py.File(path,'w') as f:
        f['costs']=np.array(rows);f.attrs['columns']='chirps,accelerations,correction_groups,bank_nodes,estimated_real_operations_per_window,required_operations_per_second,sum_requested_columns'
        f.attrs['assumptions']='4 RX; one new window per 25.37 us; complex FFT ~5 L log2 L real arithmetic operations; correction multiply 6; coarse scoring 16/node; full grid column union upper bound. Excludes phase generation, joint refinement, memory/transfer and integer indexing. Not runtime benchmark.'
    return rows
if __name__=='__main__':
    for row in costs():print(row)
