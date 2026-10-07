"""Original stored receiver I/Q for inspecting interference and clipping.

Base conda Python. No clutter subtraction, whitening, masking, conjugation,
DC removal or taper is applied to plotted raw samples or their FFTs.
The stored data are ADC counts, not calibrated analog volts.
"""
from pathlib import Path
import hashlib
import json

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from analyze_captures import DATA
from apply_measured_filter import rail_mask

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent/'output'


def metric(z):
    z = z.astype(np.complex128)
    energy = np.mean(abs(z)**2)
    repeat = np.mean(abs(z.mean(axis=0))**2)/energy
    return dict(raw_complex_rms_counts=float(np.sqrt(energy)),
                rail_valued_complex_fraction=float(rail_mask(z).mean()),
                repeating_waveform_power_fraction=float(repeat))


def decorate_wave(a, samples, fs):
    time = np.arange(len(samples))/fs*1e6
    a.plot(time, samples.real, '.-', ms=2, lw=.8, label='I')
    a.plot(time, samples.imag, '.-', ms=2, lw=.8, label='Q')
    a.axhline(-32768, color='tab:red', ls='--', lw=.6)
    a.axhline(32767, color='tab:red', ls='--', lw=.6)
    a.set(xlabel='ADC time within chirp (µs)', ylabel='Original ADC counts',
          ylim=(-35500,35500), xlim=(0,len(samples)/fs*1e6))


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.1), layout='constrained')
    overview, oa = plt.subplots(2, 2, figsize=(11,7.3), layout='constrained')
    results = {}
    with h5py.File(ROOT/'memo_004_measured_filter.h5') as prior, \
         h5py.File(OUT/'raw_receiver_iq.h5','w') as out:
        out.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        out.attrs['selection_source_sha256'] = hashlib.sha256((ROOT/'memo_004_measured_filter.h5').read_bytes()).hexdigest()
        out.attrs['processing'] = 'Original stored I/Q; no subtraction, whitening, mask, conjugation, mean removal or taper'
        out.attrs['units'] = 'ADC counts; analog voltage conversion unavailable'
        out.attrs['time_note'] = 'Within-chirp clock from recorded fs; frame-local clock assumes reader-derived 13.71us chirp interval'
        for row, name in enumerate(['parallel','perpendicular']):
            old = prior[name]; p = json.loads(old.attrs['profile_json'])
            fr = int(old.attrs['frame_of_interest'])
            report = json.loads(old[f'frame_{fr}'].attrs['report_json'])
            onset = int(report['onset_chirp'])
            starts = {'before': onset-16, 'after': onset+16}
            with h5py.File(DATA/p['file']) as source:
                frame = source['radar_cube'][0,fr]
                before_frame = source['radar_cube'][0,fr-1]
                after_frame = source['radar_cube'][0,fr+1]
            digest = hashlib.sha256(frame.tobytes()).hexdigest()
            assert digest == old[f'frame_{fr}'].attrs['frame_sha256']
            # Verify these are the same ORIGINAL samples used in the existing
            # onset fit, before any software preprocessing.
            np.testing.assert_array_equal(frame[onset:onset+8],old[f'frame_{fr}/onset_fit/observed_counts'][:])
            blocks = {key: frame[start:start+8] for key,start in starts.items()}
            shot = '6606' if name=='parallel' else '6607'
            g=out.create_group(name)
            g.attrs['source_file']=p['file'];g.attrs['source_frame_sha256']=digest
            g.attrs['source_id']=old.attrs['source_id']
            g.attrs['profile_json']=json.dumps(p);g.attrs['frame']=fr;g.attrs['data_selected_onset_chirp']=onset
            stats={key:metric(z) for key,z in blocks.items()}
            results[name]=dict(shot=shot,frame=fr,onset_chirp=onset,selected_start_chirps=starts,blocks=stats)
            g.attrs['report_json']=json.dumps(results[name])
            frequencies = np.fft.fftshift(np.fft.fftfreq(p['samples'],1/p['fs']))
            ds=g.create_dataset('frequency_hz',data=frequencies);ds.attrs['units']='Hz'
            for col,key in enumerate(['before','after']):
                z=blocks[key]
                ds=g.create_dataset(f'{key}_iq_counts',data=z,compression='gzip');ds.attrs['units']='ADC counts'
                decorate_wave(axes[row,col],z[0,0],p['fs'])
                axes[row,col].set_title(f'Shot {shot}, RX0: {key}, chirp {starts[key]}')
                axes[row,col].legend(fontsize=9,loc='upper right')
                # Rectangular, unmodified 128-sample FFT; average power over
                # eight chirps, separately per RX. No phase summation over RX.
                power=np.fft.fftshift(np.mean(abs(np.fft.fft(z,axis=-1,norm='ortho'))**2,axis=0),axes=-1)
                ds=g.create_dataset(f'{key}_raw_spectrum_counts2',data=power);ds.attrs['units']='ADC counts squared per orthonormal FFT bin'
                axes[row,2].plot(frequencies/1e6,10*np.log10(np.maximum(power[0],1e-12)),label=f'{key}: eight chirps',lw=1)
            axes[row,2].set(title=f'Shot {shot}, RX0: raw spectrum',xlabel='Frequency in stored I/Q (MHz)',
                            ylabel='Raw FFT-bin power (dB ADC-count²)',xlim=(-6.25,6.25))
            axes[row,2].legend(fontsize=9)
            # Keep a longer unmodified onset excerpt for reinspection.
            lo=max(0,onset-32);hi=min(len(frame),onset+96)
            ds=g.create_dataset('onset_excerpt_iq_counts',data=frame[lo:hi],compression='gzip');ds.attrs['units']='ADC counts';ds.attrs['start_chirp']=lo
            for z,color,label in [(before_frame,'.55',f'Pre-event frame {fr-1}'),(frame,'tab:blue',f'Marked frame {fr}'),(after_frame,'tab:orange',f'Next frame {fr+1}')]:
                rms=np.sqrt(np.mean(abs(z[:,0])**2,axis=-1))
                oa[row,0].plot(np.arange(len(rms))*p['period']*1e3,rms,color=color,label=label,lw=.8)
                ds=g.create_dataset(f'frame_{label.split()[-1]}_rx0_rms',data=rms);ds.attrs['units']='ADC counts'
            oa[row,0].axvline(onset*p['period']*1e3,color='k',ls='--',lw=.8)
            oa[row,0].set(title=f'Shot {shot}: RX0 raw RMS',xlabel='Time within each frame (ms)',ylabel='RMS complex ADC counts',ylim=(0,35000))
            oa[row,0].legend(fontsize=8)
            # Display the marked-frame raw spectrum with no software cleanup.
            sp=np.fft.fftshift(np.mean(abs(np.fft.fft(frame,axis=-1,norm='ortho'))**2,axis=1),axes=-1)
            ds=g.create_dataset('marked_frame_raw_spectrogram',data=sp,compression='gzip');ds.attrs['units']='ADC counts squared per orthonormal FFT bin; mean over four RX powers'
            im=oa[row,1].imshow(10*np.log10(np.maximum(sp.T,1e-12)),origin='lower',aspect='auto',
                extent=[0,len(frame)*p['period']*1e3,frequencies[0]/1e6,(frequencies[-1]+p['fs']/p['samples'])/1e6],
                vmin=30,vmax=100,cmap='magma')
            oa[row,1].axvline(onset*p['period']*1e3,color='cyan',ls='--',lw=.8)
            oa[row,1].set(title=f'Shot {shot}: original marked-frame spectrum',xlabel='Time within marked frame (ms)',ylabel='Stored I/Q frequency (MHz)')
            oa[row,1].grid(False);overview.colorbar(im,ax=oa[row,1],label='dB ADC-count² / FFT bin')
            # Extra plots show every RX on the same count scale.
            rxfig,rxax=plt.subplots(4,2,figsize=(10,9),layout='constrained')
            for rx in range(4):
                for col,key in enumerate(['before','after']):
                    decorate_wave(rxax[rx,col],blocks[key][0,rx],p['fs'])
                    rxax[rx,col].set_title(f'RX{rx}, {key}: chirp {starts[key]}')
                    if rx==0:rxax[rx,col].legend(fontsize=8)
            rxfig.suptitle(f'Shot {shot}: all four RX, original I/Q counts including rail-valued samples',fontsize=12)
            rxfig.savefig(OUT/f'raw_receiver_iq_{shot}_all_rx.png',dpi=160);plt.close(rxfig)
            print(name,json.dumps(results[name]),flush=True)
        out.attrs['reports_json']=json.dumps(results)
    fig.suptitle('Original stored receiver I/Q: ADC counts, including rail-valued samples\n'
                 'Before/after the data-selected change. No clutter subtraction, whitening, masking or I/Q sign change.',fontsize=12)
    fig.savefig(OUT/'raw_receiver_iq.png',dpi=180);plt.close(fig)
    overview.suptitle('Raw measured data across frames: large background already present before the event\n'
                       'Dashed lines mark the earlier data-selected onset; frame time axes are separate, with assumed chirp interval.',fontsize=11)
    overview.savefig(OUT/'raw_receiver_iq_overview.png',dpi=180);plt.close(overview)


if __name__=='__main__':
    main()
