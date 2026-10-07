"""Exclude first objectively elevated residual chirps from matched-filter search.
Uses SAME complex-row subtraction as displayed plot (0:1629 baseline).
No rescan needed: select the already exhaustively scanned time windows, then
refine the strict frame34 maximum on the same three-dimensional bank.
"""
import json,hashlib
from pathlib import Path
import h5py,numpy as np
from search_quiet_frames import QuietBank,NC,OUT
from analyze_captures import echo,C

def main():
    with h5py.File(OUT/'quiet_subtracted_6607.h5') as sub,h5py.File(OUT/'quiet_matched_filter_6607.h5') as full,h5py.File(OUT/'quiet_cutoff_audit_6607.h5','w') as out:
        z=sub['residual_complex_counts'][:1629];energy=np.mean(abs(z)**2,axis=1);train=energy[:1400]
        median=np.median(train);scale=1.4826*np.median(abs(train-median));threshold=median+6*scale
        runs=np.convolve((energy>threshold).astype(int),np.ones(3,dtype=int),'valid');onset=int(np.flatnonzero(runs>=3)[0])
        assert onset==1624
        s=full['frame_34/window_starts'][:];score=full['frame_34/coarse_score'][:];ok=s+NC<=onset
        ix=np.flatnonzero(ok)[np.argmax(score[ok])];start=int(s[ix]);p=json.loads(full.attrs['profile_json']);variance=full.attrs['reference_variance_counts2']
        bank=QuietBank(p);data=z[start:start+NC];inner=bank.q.conj()@data.ravel();inner_conj=bank.q@data.ravel();values=np.stack([abs(inner)**2,abs(inner_conj)**2])/variance
        candidates=[]
        for sign in [0,1]:
            for k in np.argsort(values[sign])[::-1][:5]:
                val,theta,fit,success=bank.refine(data.ravel(),bank.theta[k],sign,variance)
                candidates.append((val,theta,fit,success,sign))
        val,theta,fit,success,sign=max(candidates,key=lambda a:a[0]);e=echo(p,theta,NC,gate=False);doppler=-p['f_start']*2*e['radial_velocity']/(C+e['radial_velocity'])
        report=dict(frame=34,start_chirp=start,chirps=NC,coarse_score=float(score[ix]),score=val,theta=theta.tolist(),iq_conjugated=bool(sign),optimizer_success=success,explained_residual_energy_fraction=float(val*variance/np.sum(abs(data)**2)),range_m=[float(e['range'].min()),float(e['range'].max())],carrier_doppler_MHz=[float(doppler.min()/1e6),float(doppler.max()/1e6)],total_if_MHz=[float(e['if_hz'].min()/1e6),float(e['if_hz'].max()/1e6)])
        out.attrs['script']='capture_analysis/check_quiet_cutoff.py';out.attrs['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();out.attrs['report_json']=json.dumps(report);out.attrs['disturbance_first_chirp']=onset
        out.attrs['onset_rule']='First run of3chirps with mean complex residual power > median+6*1.4826*MAD of first1400chirps. Descriptive change test, not independent projectile trigger.'
        out.attrs['background_note']='Same plotted complex mean from0–1628. Audit excludes direct use of1624andlaterchirps; mean unchanged for consistency with displayed subtraction.'
        out.attrs['residual_energy_threshold_counts2']=threshold
        for k,v in [('strict_window_starts',s[ok]),('strict_coarse_score',score[ok]),('residual_chirp_energy',energy),('residual_complex_counts',data),('fitted_complex_counts',fit.reshape(NC,128)),('time_s',e['t']),('range_m',e['range']),('total_if_hz',e['if_hz']),('carrier_doppler_hz',doppler)]:out.create_dataset(k,data=v)
        print('STRICT CUTOFF',onset,'retained windows',ok.sum(),flush=True);print(json.dumps(report,indent=2),flush=True)
        print('Early full-control maximum',max(full[f'frame_{f}/coarse_score'][:].max() for f in range(34)),flush=True)
if __name__=='__main__':main()
