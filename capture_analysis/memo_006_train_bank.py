"""Automatic coherent train bank: full relative-phase bound and per-chirp FFT.
Run with base conda Python. Bounds cover x0 [m], v0 [m/s], y0 [m].
No samples are inserted in chirp gaps. No independent chirp phase is discarded.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import numpy as np
import h5py
from scipy.fft import fft, next_fast_len

C = 299792458.0


def phase(theta, t, u, f0, slope):
    x, v, y = theta
    r = np.hypot(x-v*t, y)
    tau = 2*r/C
    return 2*np.pi*(-f0*tau-slope*u*tau+.5*slope*tau*tau)


def relative_gradient_bound(box, starts, u, f0, slope):
    """Analytic Lipschitz bound for phase minus ONE train reference phase.

    Bound d(grad R)/dt and the variation of its phase multiplier separately.
    Entire parameter cells and continuous trajectory interval are covered;
    corners sampled in validation are not the source of the guarantee.
    """
    x, v, y = box
    ta, tb = starts[0]+u[0], starts[-1]+u[-1]
    tref, uref = (ta+tb)/2, (u[0]+u[-1])/2
    products = [vv*tt for vv in v for tt in [ta, tb]]
    lo, hi = x[0]-max(products), x[1]-min(products)
    xmin = 0 if lo <= 0 <= hi else min(abs(lo), abs(hi))
    ymin = 0 if y[0] <= 0 <= y[1] else min(abs(y))
    rmin = np.hypot(xmin, ymin)
    rmax = np.hypot(max(abs(lo), abs(hi)), max(abs(y)))
    if rmin == 0:
        raise ValueError('Bounds permit zero range; exclude the collision geometry')
    vm, tm = max(abs(v)), max(abs(ta), abs(tb))
    dt, du = (tb-ta)/2, (u[-1]-u[0])/2
    kmax = 4*np.pi/C*(abs(f0)+abs(slope)*max(abs(u))+2*abs(slope)*rmax/C)
    # g_x = X/R, g_v = -t X/R, g_y = y/R.
    # |d g_x/dt| <= |v|/R, |d g_y/dt| <= |v|/(2R).
    dg = dt*np.array([vm/rmin, 1+tm*vm/rmin, vm/(2*rmin)])
    dk = 4*np.pi/C*abs(slope)*(du+2*vm*dt/C)
    return (kmax*dg+dk*np.array([1, abs(tref), 1]))*(1+1e-12)


def make_bank(bounds, starts, u, f0, slope, limit, maximum):
    pending = [np.asarray(bounds, float)]
    boxes, errors = [], []
    while pending:
        box = pending.pop()
        contributions = np.diff(box, axis=1).ravel()/2*relative_gradient_bound(box, starts, u, f0, slope)
        error = contributions.sum()
        if error <= limit:
            boxes.append(box); errors.append(error)
        else:
            axis = int(np.argmax(contributions)); mid = box[axis].mean()
            if mid in box[axis]:
                raise RuntimeError('Cannot split cell at floating-point precision')
            left, right = box.copy(), box.copy()
            left[axis, 1] = right[axis, 0] = mid
            pending.extend([right, left])
        if len(boxes)+len(pending) > maximum:
            raise RuntimeError('Template limit exceeded; no incomplete bank returned')
    boxes = np.asarray(boxes)
    return boxes, boxes.mean(axis=2), np.asarray(errors)


def output(data, theta, starts, u, f0, slope, fs, nfft):
    """Complex sum from per-chirp FFTs, preserving physical frequency/phase."""
    x, v, y = theta
    uc = (u[0]+u[-1])/2; h = u-uc
    t = starts[:, None]+u
    phi = phase(theta, t, u, f0, slope)
    tc = starts+uc; r = np.hypot(x-v*tc, y)
    rd = -v*(x-v*tc)/r
    phic = phase(theta, tc, uc, f0, slope)
    frequency = -2/C*((f0+slope*uc)*rd+slope*r)+4*slope*r*rd/C**2
    correction = np.exp(1j*(phi-phic[:, None]-2*np.pi*frequency[:, None]*h))
    spectra = fft(data*np.conj(correction), nfft, axis=1)
    k = np.rint(frequency*nfft/fs).astype(np.int64)
    fk = k*fs/nfft
    values = spectra[np.arange(len(starts)), k % nfft]*np.exp(-2j*np.pi*fk*(u[0]-uc))
    return np.sum(np.exp(-1j*phic)*values)


def run(args):
    bounds = np.array([args.x0,args.v0,args.y0], float)
    if not np.all(np.isfinite(bounds)) or np.any(bounds[:,0]>bounds[:,1]):
        raise ValueError('Finite ordered parameter bounds required')
    if (not np.all(np.isfinite([args.max_loss,args.fs,args.adc_start,args.period,args.f0,args.slope]))
        or not 0 < args.max_loss < 1 or args.samples<2 or args.fs<=0
        or args.chirps<1 or args.adc_start<0 or args.max_templates<1
        or args.period <= args.adc_start+(args.samples-1)/args.fs):
        raise ValueError('Invalid clock, pulse count or loss')
    u = args.adc_start+np.arange(args.samples)/args.fs
    starts = np.arange(args.chirps)*args.period
    epsilon = np.arccos(np.sqrt(1-args.max_loss))
    boxes, centres, errors = make_bank(bounds, starts, u, args.f0, args.slope, epsilon/2, args.max_templates)
    print(f'Built {len(boxes)} complete train cells; searching {args.chirps} chirps', flush=True)
    H = (u[-1]-u[0])/2
    nfft = next_fast_len(max(args.samples, int(np.ceil(2*np.pi*args.fs*H/epsilon))))
    truth = bounds[:,0]+np.array([.37,.61,.43])*np.diff(bounds,axis=1).ravel()
    t = starts[:,None]+u
    data = np.exp(1j*phase(truth,t,u,args.f0,args.slope))
    scores = np.array([abs(output(data,p,starts,u,args.f0,args.slope,args.fs,nfft)/data.size)**2 for p in centres])
    inside = np.flatnonzero(np.all((truth>=boxes[:,:,0])&(truth<=boxes[:,:,1]),axis=1))[0]
    assert scores[inside] >= 1-args.max_loss-1e-9
    # Direct Eq.25 at every bank centre, independently of the FFT decomposition.
    direct = np.array([abs(np.sum(data*np.exp(-1j*phase(p,t,u,args.f0,args.slope)))/data.size)**2 for p in centres])
    assert direct[inside] >= np.cos(errors[inside])**2-1e-9
    # Each FFT's bin residual is bounded uniformly in fast time and therefore
    # bounds the difference in normalized COMPLEX train outputs as well.
    for i in np.linspace(0,len(centres)-1,min(32,len(centres)),dtype=int):
        z=output(data,centres[i],starts,u,args.f0,args.slope,args.fs,nfft)/data.size
        zd=np.sum(data*np.exp(-1j*phase(centres[i],t,u,args.f0,args.slope)))/data.size
        assert abs(z-zd)<=epsilon/2+1e-9
    # Check full phase, including relative chirp centres, at all 8 cell corners
    # and two random interior points of selected cells (deterministic seed).
    rng = np.random.default_rng(20261009)
    tref=(t.min()+t.max())/2; uref=(u[0]+u[-1])/2
    checked=0; largest=0.
    for i in np.linspace(0,len(boxes)-1,min(64,len(boxes)),dtype=int):
        box=boxes[i]; centre=centres[i]
        base=phase(centre,t,u,args.f0,args.slope)-phase(centre,tref,uref,args.f0,args.slope)
        points=[box[np.arange(3),corner] for corner in itertools.product([0,1],repeat=3)]
        points += [box[:,0]+rng.random(3)*np.diff(box,axis=1).ravel() for _ in range(2)]
        for p in points:
            delta=phase(p,t,u,args.f0,args.slope)-phase(p,tref,uref,args.f0,args.slope)-base
            actual=np.max(abs(delta)); assert actual<=errors[i]+1e-8
            largest=max(largest,actual);checked+=1
    volume=np.prod(np.diff(bounds,axis=1))
    np.testing.assert_allclose(np.prod(np.diff(boxes,axis=2).squeeze(-1),axis=1).sum(),volume,rtol=1e-12)
    with h5py.File(args.output,'w') as h:
        h.attrs['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        h.attrs['parameters_json']=json.dumps(vars(args),default=str)
        for k,value in [('bounds',bounds),('cells',boxes),('centres',centres),('phase_bounds',errors),('fft_retention',scores),('direct_retention',direct),('truth',truth),('chirp_start_s',starts),('fast_time_s',u),('synthetic_iq',data)]:h.create_dataset(k,data=value)
        h.attrs['nfft']=nfft;h.attrs['containing_cell']=inside;h.attrs['checked_points']=checked;h.attrs['largest_observed_cell_error']=largest
    print(f'{args.chirps} chirps: {len(boxes)} certified cells, FFT {nfft}, containing-cell power {scores[inside]:.9f}, peak {scores.max():.9f}, {checked} phase checks',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--x0',type=float,nargs=2,default=[.099,.101]);p.add_argument('--v0',type=float,nargs=2,default=[340,343]);p.add_argument('--y0',type=float,nargs=2,default=[.149,.151])
    p.add_argument('--chirps',type=int,default=4);p.add_argument('--samples',type=int,default=225);p.add_argument('--fs',type=float,default=12.5e6);p.add_argument('--adc-start',type=float,default=2e-6);p.add_argument('--period',type=float,default=25.37e-6);p.add_argument('--f0',type=float,default=77e9);p.add_argument('--slope',type=float,default=99.987e12);p.add_argument('--max-loss',type=float,default=.05);p.add_argument('--max-templates',type=int,default=200000);p.add_argument('--output',type=Path,default=Path(__file__).with_suffix('.h5'))
    run(p.parse_args())
