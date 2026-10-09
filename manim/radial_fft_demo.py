"""Numerical Doppler/acceleration illustration and direct-vs-FFT validation.

The teaching example isolates monostatic carrier phase, with uniform synthetic
sampling. The GUI uses the full sampled FMCW phase and correction groups.
Scientific numerical products are HDF5; plot PNGs are generated under media/.
"""
from pathlib import Path
import time
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits

BLUE, ORANGE, GRAY, PURPLE = '#176BB0','#B55A00','#556575','#7939A8'

def generate_demo():
    root=Path(__file__).parent
    out=root/'media/radial_fft_demo';out.mkdir(parents=True,exist_ok=True)
    n=1024;duration=144e-6;rate=n/duration;length=8*n
    wavelength=299792458./77e9;velocity=320.;acceleration=1e6
    h=(np.arange(n)-(n-1)/2)/rate
    phase_a=-2*np.pi*acceleration*h*h/wavelength
    phase_v=-4*np.pi*velocity*h/wavelength
    voltage=np.exp(1j*(phase_v+phase_a))
    frequencies=np.fft.fftfreq(length,1/rate)
    velocities=-wavelength*frequencies/2
    indices=np.flatnonzero((velocities>=0)&(velocities<=600))
    indices=indices[np.argsort(velocities[indices])];axis=velocities[indices]
    trials=np.array([0.,5e5,1e6])
    corrections=np.exp(1j*2*np.pi*trials[:,None]*h*h/wavelength)
    corrected=corrections*voltage
    coefficients=np.fft.fft(corrected,n=length,axis=1)[:,indices]
    # Independent direct inner products, including the centered-time phase.
    templates=np.exp(1j*4*np.pi*axis[:,None]*h[None,:]/wavelength)
    with threadpool_limits(limits=1):
        direct=corrected@templates.T
    shift=np.exp(1j*4*np.pi*axis*h[0]/wavelength)
    np.testing.assert_allclose(coefficients*shift,direct,rtol=1e-10,atol=2e-10)
    powers=abs(coefficients/n)**2;direct_powers=abs(direct/n)**2
    np.testing.assert_allclose(powers,direct_powers,rtol=1e-10,atol=2e-13)
    np.testing.assert_allclose(corrected[-1],np.exp(1j*phase_v),atol=2e-13)
    assert abs(axis[np.argmax(powers[-1])]-velocity)<=abs(axis[1]-axis[0])/2
    err=np.max(abs(powers-direct_powers))

    def fft_search():return abs(np.fft.fft(corrected[-1],n=length)[indices]/n)**2
    def direct_search():return abs(templates@corrected[-1]/n)**2
    with threadpool_limits(limits=1):
        for _ in range(8):fft_search();direct_search()
        fft_times=[];direct_times=[]
        for _ in range(101):
            start=time.perf_counter();direct_search();direct_times.append(time.perf_counter()-start)
            start=time.perf_counter();fft_search();fft_times.append(time.perf_counter()-start)
    fft_ms=np.median(fft_times)*1000;direct_ms=np.median(direct_times)*1000
    summary=dict(samples=n,fft_length=length,velocity_bins=len(axis),duration_s=duration,
                 synthetic_rate_hz=rate,wavelength_m=wavelength,velocity_m_s=velocity,
                 acceleration_m_s2=acceleration,maximum_power_error=err,
                 cpu_fft_ms=fft_ms,cpu_direct_ms=direct_ms,speedup=direct_ms/fft_ms)
    path=root/'assets/radial_fft_demo.h5'
    with h5py.File(path,'w') as f:
        f.attrs.update(summary);f.attrs['generator']='manim/radial_fft_demo.py'
        f.attrs['model']='synthetic uniform-sampling carrier-phase illustration, not full FMCW'
        f.attrs['benchmark']='NumPy, single BLAS thread, 101 medians; direct phasors precomputed; common correction excluded'
        for name,data in [('time_s',h),('acceleration_phase_rad',phase_a),('voltage',voltage),('corrected_voltage',corrected),('trial_acceleration_m_s2',trials),('velocity_m_s',axis),('fft_matched_power',powers),('direct_matched_power',direct_powers)]:f.create_dataset(name,data=data)

    plt.rcParams.update({'font.size':17,'axes.labelsize':18,'axes.titlesize':18,'legend.fontsize':15,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18,'savefig.facecolor':'white'})
    t=h*1e6
    def finish(fig,name):
        fig.savefig(out/(name+'.png'),dpi=180,bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(12.6,3.0),layout='constrained')
    ax.plot(t,np.cos(phase_a),color=BLUE,label='Real part')
    ax.plot(t,np.sin(phase_a),color=ORANGE,label='Imaginary part')
    ax.set(xlabel='Time from midpoint (microseconds)',ylabel='Acceleration factor',ylim=(-1.15,1.15))
    ax.legend(loc='upper center',ncol=2);finish(fig,'acceleration_factor')

    fig,axs=plt.subplots(2,1,figsize=(12.6,4.6),sharex=True,layout='constrained')
    for ax,z,title in zip(axs,[voltage,corrected[-1]],['Before correction: frequency changes with time','After the correct acceleration correction: constant frequency']):
        ax.plot(t,z.real,color=BLUE,lw=1.2,label='Real part');ax.plot(t,z.imag,color=ORANGE,lw=1.2,label='Imaginary part')
        ax.set(ylabel='Complex voltage\n(relative units)',ylim=(-1.15,1.15),title=title)
    axs[0].legend(loc='upper right',ncol=2);axs[-1].set_xlabel('Time from midpoint (microseconds)');finish(fig,'complex_voltage')

    fig,ax=plt.subplots(figsize=(12.6,4),layout='constrained')
    for a,p,col in zip(trials,powers,[GRAY,ORANGE,BLUE]):ax.plot(axis,p,color=col,lw=2,label=f'Trial acceleration = {a/1e6:g} million m/s²')
    peak=np.argmax(powers[-1]);ax.plot(axis[peak],powers[-1,peak],'o',color=BLUE)
    ax.annotate(f'Peak near {velocity:g} m/s',xy=(axis[peak],powers[-1,peak]),xytext=(405,.83),arrowprops=dict(arrowstyle='->',color=BLUE),color=BLUE)
    ax.set(xlabel='Trial radial velocity (m/s)',ylabel='Matched power (normalized)',xlim=(0,600),ylim=(0,1.12));ax.legend(loc='upper left');finish(fig,'velocity_search')

    fig,ax=plt.subplots(figsize=(7.3,4.2),layout='constrained')
    ax.plot(axis,powers[-1],color=BLUE,lw=2,label='FFT: all bins together')
    pick=np.unique(np.r_[np.arange(0,len(axis),12),np.argmax(powers[-1])])
    ax.plot(axis[pick],direct_powers[-1,pick],'o',mfc='none',color=ORANGE,ms=6,label='Direct sums: one velocity at a time')
    ax.set(xlabel='Trial radial velocity (m/s)',ylabel='Matched power (normalized)',xlim=(0,600),ylim=(0,1.08));ax.legend(loc='upper left',fontsize=12);finish(fig,'fft_equals_direct')
    # Independent conducting-sphere Mie series and all inverse branches.
    from scipy.special import spherical_jn, spherical_yn
    from scipy.optimize import brentq
    def mie(d):
        k=2*np.pi*77e9/299792458.;x=k*d/2
        order=np.arange(1,int(np.ceil(x+4*np.cbrt(x)+12))+1)
        j,y=spherical_jn(order,x),spherical_yn(order,x)
        dp=j+x*spherical_jn(order,x,derivative=True)
        dx=j+1j*y+x*(spherical_jn(order,x,derivative=True)+1j*spherical_yn(order,x,derivative=True))
        return np.pi/k**2*abs(np.sum((2*order+1)*(-1.)**order*(-dp/dx+x*j/(x*(j+1j*y)))))**2
    target=mie(.002)
    diameters=np.geomspace(.00001,.02,4097)
    sigmas=np.array([mie(d) for d in diameters])
    roots=np.array([brentq(lambda d:mie(d)-target,a,b,xtol=1e-15) for a,b,fa,fb in zip(diameters[:-1],diameters[1:],sigmas[:-1]-target,sigmas[1:]-target) if fa*fb<0])
    assert np.min(abs(roots-.002))<1e-10 and len(roots)>1
    with h5py.File(path,'a') as f:
        f['mie_diameter_m']=diameters;f['mie_rcs_m2']=sigmas;f['mie_inverse_roots_m']=roots
        f.attrs['mie_target_rcs_m2']=target
    fig,ax=plt.subplots(figsize=(12.6,4.1),layout='constrained')
    ax.plot(diameters*1000,sigmas*1e6,color=BLUE,label='Mie prediction: ideal conducting sphere, 77 GHz')
    ax.axhline(target*1e6,color=ORANGE,label='Inferred RCS (synthetic example)')
    ax.plot(roots*1000,np.full(len(roots),target*1e6),'o',color=PURPLE)
    for q in roots:ax.annotate(f'{q*1000:.3f} mm',(q*1000,target*1e6),xytext=(0,16 if q==roots[0] else -28),textcoords='offset points',ha='center',fontsize=13,color=PURPLE)
    ax.set(xlabel='Sphere diameter (mm)',ylabel='Radar cross section (mm²)',xlim=(0,max(roots)*1000+.7),ylim=(0,target*1e6*2))
    ax.legend(loc='upper left',fontsize=13);finish(fig,'mie_diameter')
    print('PASS: acceleration removed; FFT equals all direct complex sums and powers; peak within half a bin')
    print(summary)
    return summary

if __name__=='__main__':generate_demo()
