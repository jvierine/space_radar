"""Plot actual pre-disturbance Eq25 finite-bank output, not a detection claim."""
from pathlib import Path
import json,hashlib
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_captures import C,echo
from search_quiet_frames import QuietBank,NC
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'output'

def main():
    with h5py.File(OUT/'quiet_matched_filter_6607.h5') as h,h5py.File(OUT/'quiet_cutoff_audit_6607.h5') as audit,h5py.File(OUT/'quiet_match_diagnostic_6607.h5','w') as out:
        report=json.loads(audit.attrs['report_json']);p=json.loads(h.attrs['profile_json']);variance=h.attrs['reference_variance_counts2']
        cutoff=int(audit.attrs['disturbance_first_chirp'])
        z=audit['residual_complex_counts'][:];fit=audit['fitted_complex_counts'][:]
        bank=QuietBank(p);sign=report['iq_conjugated'];inner=(bank.q.conj() if not sign else bank.q)@z.ravel();surface=(abs(inner)**2/variance).reshape(*(len(g) for g in bank.grids))
        e=echo(p,report['theta'],NC,gate=False)
        doppler=-p['f_start']*2*e['radial_velocity']/(C+e['radial_velocity'])
        diagnostic={**report,'range_m':[float(e['range'].min()),float(e['range'].max())],'carrier_doppler_MHz':[float(doppler.min()/1e6),float(doppler.max()/1e6)],'total_if_MHz':[float(e['if_hz'].min()/1e6),float(e['if_hz'].max()/1e6)]}
        out.attrs['script']='capture_analysis/show_quiet_match.py';out.attrs['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();out.attrs['input_sha256']=hashlib.sha256((OUT/'quiet_matched_filter_6607.h5').read_bytes()).hexdigest();out.attrs['report_json']=json.dumps(diagnostic)
        out.attrs['strict_cutoff_chirp']=cutoff
        out.attrs['cutoff_audit_sha256']=hashlib.sha256((OUT/'quiet_cutoff_audit_6607.h5').read_bytes()).hexdigest()
        out.create_dataset('coarse_score_3d',data=surface)
        for name,g in zip(['x0_grid_m','v0_grid_m_s','y0_grid_m'],bank.grids):out.create_dataset(name,data=g)
        for k,v in [('range_m',e['range']),('radial_velocity_m_s',e['radial_velocity']),('carrier_doppler_hz',doppler),('total_if_hz',e['if_hz']),('residual_complex_counts',z),('fitted_complex_counts',fit)]:out.create_dataset(k,data=v)
        frames=list(range(35));maxima=[h[f'frame_{f}/coarse_score'][:].max() for f in frames]
        contaminated_max=maxima[-1];maxima[-1]=audit['strict_coarse_score'][:].max()
        fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
        a=ax[0,0];a.plot(frames,maxima,'o-',ms=3,label='Strict quiet windows');a.plot(34,contaminated_max,'rx',ms=8,label='Includes first disturbed chirps');a.legend(fontsize=9);a.set(xlabel='Frame index',ylabel='Maximum matched-filter score',title='All pre-disturbance windows; both I/Q signs')
        a=ax[0,1];s=h['frame_34/window_starts'][:];vals=h['frame_34/coarse_score'][:];a.plot(s,vals,lw=.8);a.axvspan(cutoff-NC+1,1622,color='r',alpha=.18,label='Excluded: includes chirp 1624 or later');a.plot(report['start_chirp'],report['coarse_score'],'o',color='orange');a.legend(fontsize=8);a.set(xlabel='First chirp of 8-chirp window in frame 34',ylabel='Maximum score over x0, v0, y0',title='Strict windows end before chirp 1624')
        im=ax[1,0].pcolormesh(bank.grids[1]/1000,bank.grids[0],surface.max(axis=2),shading='auto',cmap='viridis');ax[1,0].set(xlabel='v0 (km/s)',ylabel='x0 (m)',title=f'Best window: frame {report["frame"]}, chirp {report["start_chirp"]}; max over y0');fig.colorbar(im,ax=ax[1,0],label='Matched-filter score')
        y=surface.max(axis=(0,1));ax[1,1].plot(bank.grids[2],y,'o-');ax[1,1].set(xlabel='y0 (m)',ylabel='Maximum score over x0, v0',title='Profiled transverse-distance scan')
        fig.suptitle('Shot 6607 — matched filter before sustained disturbance at chirp 1624\nScore is a residual-variance-normalized projection, not calibrated thermal SNR',fontsize=13)
        fig.text(.01,.001,'Script: capture_analysis/show_quiet_match.py; bank: search_quiet_frames.py',fontsize=7,color='gray');fig.savefig(OUT/'quiet_matched_filter_6607.png',dpi=180);plt.close(fig)
        fig,ax=plt.subplots(2,1,figsize=(12,6.5),layout='constrained');t=e['t']*1e6
        for k in range(NC):
            ax[0].plot(t[k],z[k].real,color='gray',lw=.7,label='Real residual' if k==0 else None);ax[0].plot(t[k],fit[k].real,color='tab:red',lw=.8,label='Real fitted template' if k==0 else None)
            ax[1].plot(t[k],e['if_hz'][k]/1e6,color='k');ax[1].plot(t[k],doppler[k]/1e6,color='tab:blue')
        ax[0].legend();ax[0].set(ylabel='ADC counts',title=f'Frame {report["frame"]}, chirps {report["start_chirp"]}–{report["start_chirp"]+7}: fit explains {report["explained_residual_energy_fraction"]*100:.2f}% of residual energy')
        ax[1].plot([],[],color='k',label='Total physical IF');ax[1].plot([],[],color='tab:blue',label='Carrier Doppler');ax[1].legend();ax[1].set(xlabel='Time within fitted window (µs)',ylabel='Signed frequency (MHz)')
        fig.suptitle(f'Conditional fit: x0={report["theta"][0]:.4f} m, v0={report["theta"][1]/1000:.4f} km/s, y0={report["theta"][2]:.4f} m\nNo confirmed pellet echo; unknown ADC sign/sideband and chamber geometry',fontsize=12)
        fig.text(.01,.001,'Script: capture_analysis/show_quiet_match.py',fontsize=7,color='gray');fig.savefig(OUT/'quiet_match_waveform_6607.png',dpi=180);plt.close(fig)
        print(json.dumps(diagnostic,indent=2),flush=True)
if __name__=='__main__':main()
