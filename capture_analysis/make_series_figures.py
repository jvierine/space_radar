"""Figures and scoped HDF5 inputs for the four chamber radar memos.

All figures use a final width of 7.06 inches and >=10 pt labels. Numerical
inputs, assumptions and script hashes stay in HDF5, never CSV.
"""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_captures import C, KB, TSYS, PT, G, echo, sphere

ROOT = Path(__file__).resolve().parent
FIG = ROOT.parent/'figures'
NAMES = ['parallel','perpendicular']
LABELS = ['6606: parallel profile','6607: perpendicular profile']


def provenance(h):
    h.attrs['figure_script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    h.attrs['model_script_sha256'] = hashlib.sha256((ROOT/'analyze_captures.py').read_bytes()).hexdigest()
    h.attrs['simulation_status'] = 'Synthetic complex I/Q; recorded profiles, assumed geometry/material/gains/filtering'


def dataset(g,name,value,units):
    ds=g.create_dataset(name,data=value,compression='gzip' if np.ndim(value) else None)
    ds.attrs['units']=units


def gap(t, z):
    return (np.pad(t,((0,0),(0,1)),constant_values=np.nan).ravel(),
            np.pad(z,((0,0),(0,1)),constant_values=np.nan).ravel())


def save(fig,name):
    fig.savefig(FIG/name)
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.titlesize':10,'axes.labelsize':10,
        'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,
        'axes.grid':True,'grid.alpha':.2})
    with h5py.File(ROOT/'search_experiment.h5') as sim, h5py.File(ROOT/'feasibility.h5') as old, \
         h5py.File(ROOT/'memo_001_signal.h5','w') as signal, h5py.File(ROOT/'memo_002_snr.h5','w') as snr:
        provenance(signal); provenance(snr)
        signal.attrs['source_simulation_sha256']=hashlib.sha256((ROOT/'search_experiment.h5').read_bytes()).hexdigest()
        signal.attrs['assumptions_json']=old.attrs['assumptions_json']
        snr.attrs['assumptions_json']=old.attrs['assumptions_json']
        waves, axes=plt.subplots(2,2,figsize=(7.06,4.7),layout='constrained')
        freq, fax=plt.subplots(1,2,figsize=(7.06,2.9),layout='constrained')
        sf, sax=plt.subplots(1,2,figsize=(7.06,3.0),layout='constrained')
        counts=np.array([1,2,4,8,16,32])
        for row,(name,label) in enumerate(zip(NAMES,LABELS)):
            g=sim[name]; p=json.loads(g.attrs['profile_json']); theta=g['truth'][:]
            e=echo(p,theta,32)
            exact=g['exact_echo'][:]; noisy=g['noisy_iq'][:]
            np.testing.assert_allclose(exact,e['iq'],rtol=1e-12,atol=1e-18)
            sg=signal.create_group(name); sg.attrs['profile_json']=json.dumps(p)
            dataset(sg,'theta',theta,'columns: m, m/s, m')
            for key,units in [('t','s'),('u','s'),('range','m'),('radial_velocity','m/s'),
                              ('if_hz','Hz'),('phase','rad'),('active','1'),('retained','1')]:
                dataset(sg,key,e[key],units)
            dataset(sg,'exact_iq',exact,'sqrt(W), assumed input-equivalent power-wave scale')
            dataset(sg,'noisy_iq',noisy,'sqrt(W), same assumed scale')
            dataset(sg,'noise_variance',g['noise_variance_W'][()],'W')
            first=0 if name=='parallel' else 6
            sg.attrs['matched_filter_window_start_chirp']=first
            norm=max(abs(exact).max(),1e-30)
            ta,za=gap(e['t']*1e6,abs(exact)/norm)
            axes[row,0].plot(ta,za,color='k')
            axes[row,0].set(xlabel='Time from path start (µs)',ylabel='Normalized |s|',title=label)
            axes[row,0].set_xlim(0,175)
            idx=np.arange(p['samples'])/p['fs']<2e-6
            u=np.arange(p['samples'])[idx]/p['fs']*1e6
            axes[row,1].plot(u,exact[first,idx].real/norm,'.-',label='I')
            axes[row,1].plot(u,exact[first,idx].imag/norm,'.-',label='Q')
            axes[row,1].set(xlabel='ADC time in chirp (µs)',ylabel='Normalized I/Q',title=f'Noiseless samples, chirp {first}')
            axes[row,1].legend(loc='upper right')
            active=e['active']
            fax[row].axhspan(-10,0,color='tab:green',alpha=.12)
            fax[row].scatter(e['t'][active]*1e6,e['if_hz'][active]/1e6,s=3,rasterized=True)
            fax[row].axhline(0,color='k',lw=.5)
            fax[row].set(xlabel='Time from path start (µs)',ylabel='Physical IF (MHz)',title=label,xlim=(0,175))
            ss=snr.create_group(name); ss.attrs['profile_json']=json.dumps(p)
            sigma=float(abs(sphere(p['f_start'],p['projectile_size']*1e-3))**2)
            pr=float(PT*G*G*(C/p['f_start'])**2*sigma/(4*np.pi)**3)
            reference=pr*counts*p['samples']/p['fs']/(KB*TSYS)
            en=np.sum(abs(exact)**2,axis=1)
            best=np.array([np.max(np.convolve(en,np.ones(n),'valid'))/p['fs']/(KB*TSYS) for n in counts])
            for key,value,unit in [('rcs',sigma,'m^2'),('power_1m',pr,'W'),('chirps',counts,'1'),
                ('adc_on_time',counts*p['samples']/p['fs'],'s'),
                ('wall_span',(counts-1)*p['period']+p['samples']/p['fs'],'s'),
                ('constant_range_snr',reference,'1'),('moving_path_snr',best,'1')]:
                dataset(ss,key,value,unit)
            sax[row].plot(counts,10*np.log10(reference),'o-',label='Constant 1 m range')
            sax[row].plot(counts,10*np.log10(best),'s-',label='Moving 1 m path')
            sax[row].set(xscale='log',xticks=counts,xticklabels=[str(i) for i in counts],
                         xlabel='Integrated chirps',ylabel='Conditional thermal SNR (dB)',title=label)
            sax[row].legend(loc='lower right')
        save(waves,'memo_001_waveform.pdf'); save(freq,'memo_001_if.pdf'); save(sf,'memo_002_snr.pdf')
        # These are precisely the synthetic records that memo 3 filters.
        g=sim['perpendicular']; p=json.loads(g.attrs['profile_json']); noise=np.sqrt(g['noise_variance_W'][()])
        fig,ax=plt.subplots(4,2,figsize=(7.06,6.6),layout='constrained')
        for row,n in enumerate([1,2,4,8]):
            case=g[f'case_{n}']; report=json.loads(case.attrs['report_json']); start=report['first_chirp']
            t=np.arange(n)[:,None]*p['period']+np.arange(p['samples'])[None,:]/p['fs']
            for col,component in enumerate([np.real,np.imag]):
                for z,l,c,w in [(g['noisy_iq'][start:start+n],'Noisy','0.65',.65),
                                (g['exact_echo'][start:start+n],'Injected','k',.65),
                                (case['fitted_echo'][:],'Fitted','tab:blue',1.)]:
                    tt,zz=gap(t*1e6,component(z)/noise)
                    ax[row,col].plot(tt,zz,label=l,color=c,lw=w)
                ax[row,col].set(xlabel='Window time (µs)',ylabel=('I' if col==0 else 'Q')+' / σ',title=f'{n} chirp'+('s' if n>1 else ''))
                if row==0 and col==0: ax[row,col].legend(loc='upper right')
        save(fig,'memo_003_recovery.pdf')
        fig,ax=plt.subplots(1,2,figsize=(7.06,3.1),layout='constrained')
        for a,name,label in zip(ax,NAMES,LABELS):
            g=sim[name]; grid=g['case_8/coarse_scores'][:].max(axis=2)
            im=a.imshow(grid,origin='lower',aspect='auto',extent=[g['v0_grid'][0],g['v0_grid'][-1],g['x0_grid'][0],g['x0_grid'][-1]])
            a.grid(False); a.set(xlabel='v₀ (m/s)',ylabel='x₀ (m)',title=label)
            fig.colorbar(im,ax=a,label='Fixed-bank score',pad=.02)
        save(fig,'memo_003_search.pdf')
    with h5py.File(ROOT/'memo_004_measured_filter.h5') as measured:
        fig,ax=plt.subplots(2,2,figsize=(7.06,5.0),layout='constrained')
        fitfig,fitax=plt.subplots(2,2,figsize=(7.06,4.5),layout='constrained')
        for row,(name,label) in enumerate(zip(NAMES,LABELS)):
            g=measured[name]; p=json.loads(g.attrs['profile_json']); reports=json.loads(g.attrs['reports_json'])
            event=next(r for r in reports if r['role']=='event'); roi=event['frame']
            for fr,c,l in [(roi-1,'0.55','Before'),(roi,'tab:blue','Marked frame'),(roi+1,'tab:orange','After')]:
                en=g[f'frame_{fr}/band_energy'][:]
                ax[row,0].plot(np.arange(len(en))*p['period']*1e3,en,color=c,label=f'{l}: {fr}',lw=.7)
            ax[row,0].axvline(event['onset_chirp']*p['period']*1e3,color='k',ls='--',lw=.8)
            ax[row,0].axvline(event['start_chirp']*p['period']*1e3,color='tab:red',ls=':',lw=.8)
            ax[row,0].set(yscale='log',xlabel='Time within each frame (ms)',ylabel='Whitened band energy',title=label)
            ax[row,0].legend(loc='upper left',fontsize=10)
            spec=g[f'frame_{roi}/whitened_spectrogram'][:]
            im=ax[row,1].imshow(10*np.log10(np.maximum(np.fft.fftshift(spec,axes=1).T,1e-12)),origin='lower',aspect='auto',
                extent=[0,len(spec)*p['period']*1e3,-p['fs']/2e6,p['fs']/2e6],vmin=0,vmax=30,cmap='magma')
            ax[row,1].grid(False)
            ax[row,1].set(xlabel='Time in marked frame (ms)',ylabel='Sampled frequency (MHz)',title=f'Measured frame {roi}')
            fig.colorbar(im,ax=ax[row,1],label='Whitened power (dB)',pad=.02)
            og=g[f'frame_{roi}/onset_fit']
            observed=og['clutter_subtracted_counts'][:]; fitted=og['fitted_counts'][:]
            good=g['valid_adc_mask'][0].astype(bool)
            t=np.arange(8)[:,None]*p['period']+np.arange(p['samples'])[None,:]/p['fs']
            # RX0, first two chirps, with invalid/clipped-training cells excluded.
            for col,component in enumerate([np.real,np.imag]):
                for z,l,c in [(observed,'Measured','0.6'),(fitted,'Fitted Eq.25','tab:blue')]:
                    values=component(z[:2,0]).copy(); values[:,~good]=np.nan
                    tt,zz=gap(t[:2]*1e6,values)
                    fitax[row,col].plot(tt,zz,color=c,label=l,lw=.8)
                fitax[row,col].set(xlabel='Time from onset window (µs)',ylabel=('I' if col==0 else 'Q')+' (ADC counts)',title=label+' — RX0')
                if row==0 and col==0: fitax[row,col].legend(loc='upper right')
        save(fig,'memo_004_diagnostics.pdf'); save(fitfig,'memo_004_onset_fit.pdf')
    print('Generated four-memo figures and scoped signal/SNR HDF5 products in',FIG)


if __name__=='__main__': main()
