"""Independent quadratic-phase reference against shipped Rust/Wasm FFT search.
Run with conda run -n base python tools/validate_radial_search.py.
Only synthetic arrays are saved, to HDF5.
"""
from pathlib import Path
import json, subprocess, time, hashlib
import numpy as np
import h5py
ROOT=Path(__file__).resolve().parents[1]
p=dict(rows=1000,n=225,frame=125,fs=12.5e6,adc=2e-6,period=25.37e-6,f0=77e9,slope=99.98738765716553e12)
C=299792458.
def template(th,n):
    u=p['adc']+np.arange(p['n'])/p['fs']
    t=np.arange(n)[:,None]*p['period']+u
    star=(t[0,0]+t[-1,-1])/2
    h=t-star;r=th[0]+th[1]*h+.5*th[2]*h*h
    tau=2*r/C
    return np.exp(2j*np.pi*(-(p['f0']+p['slope']*u)*tau+.5*p['slope']*tau*tau))*(r>0)*(u>=tau)

def run(n,bounds,truth,orientation=0,background_amplitude=0):
    rng=np.random.default_rng(638+n+orientation)
    z=(rng.normal(size=(1000,225))+1j*rng.normal(size=(1000,225)))*20/np.sqrt(2)
    z+=background_amplitude*np.exp(2j*np.pi*np.arange(225)/23)
    q=template(truth,n)
    z[750:750+n]+=q.conj() if orientation else q
    # Strong coherent injection; no use of Wasm to generate the expected waveform.
    z[750:750+n]+=5*(q.conj() if orientation else q)
    z[880:884]=0
    z=z.astype(np.complex64)
    result=subprocess.run(['node','tools/benchmark_radial.mjs'],input=json.dumps(dict(p=p,data=np.stack([z.real,z.imag],-1).ravel().tolist(),pulses=n,bounds=bounds)),capture_output=True,text=True,check=True,cwd=ROOT)
    print(result.stderr,flush=True)
    row=json.loads(result.stdout);b=np.asarray(row['best']);nr,na,nv=map(int,row['info'][:3]);scores=np.asarray(row['scores']).reshape(nr,na,nv)
    assert np.count_nonzero(np.isfinite(scores))>0
    assert all(bounds[2*i]-1e-4<=b[i]<=bounds[2*i+1]+1e-4 for i in range(3))
    assert int(b[4])==orientation,(b,orientation)
    residual=z.astype(np.complex128)-z[:125].astype(np.complex128).mean(0)
    event=residual[750:750+n].ravel();reference=template(b[:3],n).ravel()
    if orientation:reference=reference.conj()
    starts=[];k=125
    while k+n<=625:
        if k//125==(k+n-1)//125: starts.append(k);k+=1
        else:k+=1
    quiet=np.stack([z[k:k+n].ravel() for k in starts])
    # prepare chooses at most 24 approximately uniformly spaced intact controls.
    indices=np.floor(np.arange(min(24,len(quiet)))*len(quiet)/min(24,len(quiet))).astype(int)
    quiet=quiet[indices]
    exact=abs(np.vdot(reference,event))**2/np.mean(abs(quiet@reference.conj())**2)
    np.testing.assert_allclose(exact,b[3],rtol=2e-4)
    expected=reference*np.vdot(reference,event)/np.vdot(reference,reference)
    fit=np.asarray(row['fitted']).reshape(-1,2);fit=fit[:,0]+1j*fit[:,1]
    np.testing.assert_allclose(fit,expected,rtol=3e-4,atol=.003)
    coherence=abs(np.vdot(reference,q.conj().ravel() if orientation else q.ravel()))**2/(np.vdot(reference,reference).real*np.vdot(q.ravel(),q.ravel()).real)
    # Raw reference can contain narrowband stationary clutter. Its colored
    # denominator can move the ratio maximum away from the injected waveform;
    # verify the exact raw-reference statistic, not trajectory recovery, there.
    if background_amplitude==0:assert coherence>.85,(n,b,coherence)
    row['background_amplitude']=background_amplitude
    row['coherence']=coherence
    # Every valid node is searched, and both requested MAX projections retain the full bounds.
    row['projection_v_r']=np.nanmax(scores,axis=1).tolist()
    row['projection_v_a']=np.nanmax(scores,axis=0).tolist()
    return row

out=ROOT/'web/lab/qa/radial_search_validation.h5';out.parent.mkdir(parents=True,exist_ok=True)
with h5py.File(out,'w') as h:
    h.attrs['generator']='tools/validate_radial_search.py + tools/benchmark_radial.mjs'
    h.attrs['wasm_sha256']=hashlib.sha256((ROOT/'web/lab/core.wasm').read_bytes()).hexdigest()
    for name,n,bounds,truth,orient,background_amplitude in [
        ('default_8',8,[.05,.7,-400,400,0,1e6],[.3,-320,5e5],0,0),
        ('aliased_conjugate_8',8,[.28,.32,280,360,3e5,7e5],[.3,320,5e5],1,0),
        ('focused_16',16,[.28,.32,-340,-300,4e5,6e5],[.3,-320,5e5],0,0),
        ('raw_clutter_reference_8',8,[.05,.7,-400,400,0,1e6],[.3,-320,5e5],0,400),
    ]:
        row=run(n,bounds,truth,orient,background_amplitude);g=h.create_group(name)
        g.attrs['bounds_r_v_a']=bounds;g.attrs['seconds']=row['seconds'];g.attrs['coherence']=row['coherence']
        for key in ['info','best','scores','fitted','projection_v_r','projection_v_a']:g.create_dataset(key,data=np.asarray(row[key]),compression='gzip')
        print(name,'PASS, coherence',row['coherence'],flush=True)
print('PASS:',out)
