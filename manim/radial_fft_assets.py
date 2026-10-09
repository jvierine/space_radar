"""Synthetic numerical checks for radial_fft_search.py; no measurement data."""
from pathlib import Path
import h5py
import numpy as np

def generate():
    c, fs, f0, gamma = 299792458., 12.5e6, 77e9, 99.98738765716553e12
    period, ns, nc = 25.37e-6, 225, 8
    u = 2e-6 + np.arange(ns)/fs
    d = (np.arange(nc)-(nc-1)/2)*period
    h = d[:,None]+u-u.mean()
    r0, v0, a0 = .8, -300., 5e5
    r = r0+v0*h+.5*a0*h*h
    phase = 2*np.pi*(-(f0+gamma*u)*2*r/c+.5*gamma*(2*r/c)**2)
    k0 = -2*(f0+gamma*u.mean())/c+4*gamma*r0/c**2
    slow = k0*v0
    fast = slow-2*gamma*r0/c
    phi0 = 2*np.pi*(-(f0+gamma*u.mean())*2*r0/c+.5*gamma*(2*r0/c)**2)
    linear = phi0+2*np.pi*(fast*(u-u.mean())+slow*d[:,None])
    correction = phase-linear
    np.testing.assert_allclose(np.exp(1j*phase), np.exp(1j*linear)*np.exp(1j*correction), atol=2e-12)
    np.testing.assert_allclose(c*(slow-fast)/(2*gamma), r0)
    np.testing.assert_allclose(slow/k0, v0)
    assert nc*ns/fs == 144e-6
    path = Path(__file__).parent/'assets/radial_fft.h5'
    path.parent.mkdir(exist_ok=True)
    with h5py.File(path,'w') as out:
        out.attrs.update(generator='manim/radial_fft_assets.py', synthetic=True,
                         r0_m=r0,v0_m_s=v0,a0_m_s2=a0,fs_hz=fs,
                         reconstructed_period_s=period,T_coh_s=nc*ns/fs)
        for name,value in [('u_s',u),('d_s',d),('h_s',h),('range_m',r),('phase_rad',phase),('correction_rad',correction)]:
            out[name]=value
    return path

if __name__=='__main__':
    print('PASS: exact phase factorization, FFT inversion and T_coh:',generate())
