"""Signed velocity response of measured quiet frame34: no fast-echo claim."""
from pathlib import Path
import hashlib,json
import h5py,numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
from search_signed_velocities import SignedBank,OUT,NC
from analyze_captures import echo,C

def main():
    with threadpool_limits(8),torch.no_grad(),h5py.File(OUT/'quiet_signed_velocity_6607.h5') as h,h5py.File(OUT/'quiet_subtracted_6607.h5') as sub,h5py.File(OUT/'quiet_signed_velocity_diagnostic_6607.h5','w') as out:
        p=json.loads(h.attrs['profile_json']);g=h['frame_34'];report=json.loads(g.attrs['best_refined_report_json']);variance=h.attrs['reference_variance_counts2'];cutoff=int(h.attrs['strict_cutoff_chirp'])
        bank=SignedBank(p);z=sub['residual_complex_counts'][:cutoff];starts=np.arange(len(z)-NC+1);profile=np.zeros(71)
        for a in range(0,len(starts),256):
            d=np.stack([z[s:s+NC].ravel() for s in starts[a:a+256]]).T.astype(np.complex64);powers=bank.accelerated(d)
            val=torch.maximum(powers[0],powers[1]).reshape(49,71,20,-1).amax(dim=(0,2,3)).cpu().numpy()/variance;profile=np.maximum(profile,val)
        profile=np.concatenate([profile[1:][::-1],profile])
        surface=g['best_signed_coarse_surface'][:];vel=h['v0_grid_m_s'][:];x=h['x0_grid_m'][:]
        e=echo(p,report['theta'],NC,gate=False);doppler=-p['f_start']*2*e['radial_velocity']/(C+e['radial_velocity']);rangebeat=-p['freq_slope']*2*e['range']/C
        diagnostic={**report,'range_m':[float(e['range'].min()),float(e['range'].max())],'carrier_doppler_MHz':[float(doppler.min()/1e6),float(doppler.max()/1e6)],'total_if_MHz':[float(e['if_hz'].min()/1e6),float(e['if_hz'].max()/1e6)]}
        out.attrs['script']='capture_analysis/show_signed_velocity.py';out.attrs['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();out.attrs['input_sha256']=hashlib.sha256((OUT/'quiet_signed_velocity_6607.h5').read_bytes()).hexdigest();out.attrs['report_json']=json.dumps(diagnostic)
        out.attrs['profile_scope']='Maximum over all1617quiet8chirpwindows,x0,y0,bothIQsigns; signed velocities100m/sgrid,mirror symmetryexact.'
        for key,v in [('v0_grid_m_s',vel),('maximum_score_vs_velocity',profile),('best_signed_coarse_surface',surface),('range_m',e['range']),('carrier_doppler_hz',doppler),('range_beat_hz',rangebeat),('total_if_hz',e['if_hz'])]:out.create_dataset(key,data=v)
        fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
        ax[0,0].plot(vel/1000,profile);ax[0,0].axvline(0,color='gray',ls=':');ax[0,0].set(xlabel='Signed along-path v0 (km/s)',ylabel='Maximum matched-filter score',title='Quiet frame34: maximum over time, x0, y0, I/Q sign',xlim=(-7,7))
        ax[0,1].plot(g['window_starts'][:],g['coarse_score'][:],lw=.8);ax[0,1].plot(report['start_chirp'],report['coarse_score'],'o',color='orange');ax[0,1].set(xlabel='Start chirp of eight-chirp window',ylabel='Maximum score over x0, v0, y0',title='All windows end before first disturbed chirp1624')
        m=ax[1,0].pcolormesh(vel/1000,x,surface.max(2),shading='auto',cmap='viridis');ax[1,0].set(xlabel='Signed along-path v0 (km/s)',ylabel='x0 (m)',title=f'Frame34, chirp{report["start_chirp"]}: maximum over y0');fig.colorbar(m,ax=ax[1,0],label='Matched-filter score')
        frames=np.arange(35);values=[h[f'frame_{k}/coarse_score'][:].max() for k in frames];ax[1,1].plot(frames,values,'o-',ms=3);ax[1,1].set(xlabel='Frame index',ylabel='Maximum matched-filter score',title='Earlier quiet-data maxima are larger')
        fig.suptitle('Shot6607 — broad signed-velocity sweep of quiet complex residual\nFinite bank; scalar residual-variance normalization, not calibrated thermal SNR',fontsize=13)
        fig.text(.01,.001,'Script: capture_analysis/show_signed_velocity.py; bank: search_signed_velocities.py',fontsize=7,color='gray');fig.savefig(OUT/'quiet_signed_velocity_6607.png',dpi=180);plt.close(fig)
        data=g['best_residual_complex_counts'][:];fit=g['best_fitted_complex_counts'][:];t=e['t']*1e6
        fig,ax=plt.subplots(2,1,figsize=(12,7),layout='constrained')
        for k in range(NC):
            ax[0].plot(t[k],data[k].real,color='gray',lw=.8,label='Real quiet residual' if k==0 else None);ax[0].plot(t[k],fit[k].real,color='r',label='Real fitted template' if k==0 else None)
            ax[1].plot(t[k],doppler[k]/1e6,color='tab:blue',label='Carrier Doppler' if k==0 else None);ax[1].plot(t[k],rangebeat[k]/1e6,color='tab:green',label='Range beat' if k==0 else None);ax[1].plot(t[k],e['if_hz'][k]/1e6,color='k',label='Total physical IF' if k==0 else None)
        ax[0].legend();ax[0].set(ylabel='ADC counts',title=f'Frame34, chirps{report["start_chirp"]}–{report["start_chirp"]+7}: fit explains {100*report["explained_residual_energy_fraction"]:.2f}% of energy')
        ax[1].legend();ax[1].set(xlabel='Time within fitted window (µs)',ylabel='Signed frequency (MHz)')
        fig.suptitle(f'Conditional template fit: x0={report["theta"][0]:.3f}m, v0={report["theta"][1]/1000:.3f}km/s, y0={report["theta"][2]:.3f}m\nMirror(-x0,-v0) identical; parameters are not evidence of an object at that speed',fontsize=12)
        fig.text(.01,.001,'Script: capture_analysis/show_signed_velocity.py',fontsize=7,color='gray');fig.savefig(OUT/'quiet_signed_velocity_waveform_6607.png',dpi=180);plt.close(fig)
        print(json.dumps(diagnostic,indent=2),flush=True)
if __name__=='__main__':main()
