"""Separate carrier Doppler, range beat, and chirp-motion term in Eq.25.

Base conda Python. Existing measured-data fits are audited, not accepted
as fast echoes. Plotting inputs and physical sign conventions remain in HDF5.
"""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_captures import C, echo

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent/'output'


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
    fig, ax = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    with h5py.File(ROOT/'memo_004_measured_filter.h5') as source, h5py.File(OUT/'doppler_audit.h5','w') as out:
        out.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        out.attrs['input_sha256'] = hashlib.sha256((ROOT/'memo_004_measured_filter.h5').read_bytes()).hexdigest()
        out.attrs['sign_convention'] = 'rx*conj(LO); radial velocity positive receding; Doppler negative receding; stored I/Q conjugation is separate'
        for row, name in enumerate(['parallel', 'perpendicular']):
            g = source[name]; p = json.loads(g.attrs['profile_json'])
            report = json.loads(g[f"frame_{g.attrs['frame_of_interest']}"].attrs['report_json'])
            e = echo(p, report['recovered_theta'], 8)
            tau = 2*e['range']/C
            derivative = 2*e['radial_velocity']/(C+e['radial_velocity'])
            doppler = -p['f_start']*derivative
            range_beat = -p['freq_slope']*tau
            chirp_motion = -p['freq_slope']*(e['u']-tau)*derivative
            np.testing.assert_allclose(doppler+range_beat+chirp_motion, e['if_hz'], rtol=1e-12, atol=1e-7)
            mask = abs(e['iq'])>0
            scale = np.sign(e['radial_velocity'])
            label = '6606' if name == 'parallel' else '6607'
            expected = 2*p['f_start']*p['speed']/C
            a = ax[row,0]
            b = ax[row,1]
            for ch in range(8):
                t = e['t'][ch]*1e6
                a.plot(t, e['radial_velocity'][ch]/1e3, color='.6', lw=1)
                a.plot(t, np.where(mask[ch], e['radial_velocity'][ch]/1e3, np.nan), color='tab:orange', lw=2)
                for z,c in [(doppler,'tab:blue'),(range_beat,'tab:green'),(e['if_hz'],'k')]:
                    b.plot(t, z[ch]/1e6, color=c, lw=1)
                if mask[ch].any():
                    b.axvspan(t[mask[ch]].min(),t[mask[ch]].max(),color='tab:orange',alpha=.15)
            a.axhline(p['speed']/1e3,color='tab:red',ls='--',lw=.8)
            a.axhline(-p['speed']/1e3,color='tab:red',ls='--',lw=.8)
            a.plot([],[],color='.6',label='Fitted radial velocity')
            a.plot([],[],color='tab:orange',lw=2,label='Retained by previous template')
            a.plot([],[],color='tab:red',ls='--',label='Recorded projectile speed ±v')
            a.set(title=f'Shot {label}: fit passes close to radar',xlabel='Time in fitted window (µs)',ylabel='Radial velocity (km/s)',ylim=(-7.4,7.4))
            a.legend(fontsize=8,loc='upper right')
            b.axhline(expected/1e6,color='tab:red',ls='--',lw=.8)
            b.axhline(-expected/1e6,color='tab:red',ls='--',lw=.8)
            for c,text in [('tab:blue','Carrier Doppler'),('tab:green','Range beat'),('k','Total IF (includes chirp motion)')]:
                b.plot([],[],color=c,label=text)
            b.plot([],[],color='tab:red',ls='--',label='Expected ± full radial Doppler')
            b.set(title=f'Shot {label}: physical Eq.25 frequencies',xlabel='Time in fitted window (µs)',ylabel='Frequency (MHz)',ylim=(-4.1,4.1))
            b.legend(fontsize=8,loc='lower left')
            result=dict(expected_full_radial_MHz=float(expected/1e6),expected_at_10deg_MHz=float(expected*np.cos(np.deg2rad(10))/1e6),
                previous_fit_retained_carrier_doppler_MHz=[float(doppler[mask].min()/1e6),float(doppler[mask].max()/1e6)],
                previous_fit_retained_range_beat_MHz=[float(range_beat[mask].min()/1e6),float(range_beat[mask].max()/1e6)],
                retained_chirps=np.flatnonzero(mask.any(axis=1)).tolist(),iq_conjugated=report['iq_conjugated'])
            group=out.create_group(name);group.attrs['report_json']=json.dumps(result)
            for key,z,unit in [('sample_time',e['t'],'s'),('radial_velocity',e['radial_velocity'],'m/s'),
                    ('carrier_doppler',doppler,'Hz'),('range_beat',range_beat,'Hz'),('chirp_motion',chirp_motion,'Hz'),
                    ('total_if',e['if_hz'],'Hz'),('retained_mask',mask,'1')]:
                ds=group.create_dataset(key,data=z,compression='gzip');ds.attrs['units']=unit
            print(name,json.dumps(result),flush=True)
    fig.suptitle('Doppler audit: prior fits select low radial velocity near closest approach\n'
                 'Orange segments are the only samples used by those templates; physical sign is approaching +, receding −.',fontsize=11)
    fig.savefig(OUT/'doppler_audit.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    main()
