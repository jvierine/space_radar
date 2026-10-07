"""Check all recorded frames for return to the pre-event real waveform.

Reference: RX0 real waveform averaged over the first 512 chirps in the
marked frame. Compare RMS difference in nonoverlapping 128-chirp blocks.
The descriptive baseline bound is the 99th percentile over preceding frames;
it is not a calibrated hypothesis test. Original samples remain unchanged.
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

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / 'output'
BLOCK = 128


def main():
    with h5py.File(OUT / 'raw_receiver_iq.h5', 'r') as prior, \
         h5py.File(OUT / 'raw_recovery.h5', 'w') as product:
        product.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        product.attrs['method'] = __doc__
        for name, shot in [('parallel', '6606'), ('perpendicular', '6607')]:
            previous = prior[name]
            p = json.loads(previous.attrs['profile_json'])
            marked = int(previous.attrs['frame'])
            onset = int(previous.attrs['data_selected_onset_chirp'])
            rms, repeat, residual = [], [], []
            with h5py.File(DATA / p['file'], 'r') as source:
                cube = source['radar_cube']
                reference = cube[0, marked, :512, 0, :].real.astype(np.float64).mean(axis=0)
                for fr in range(cube.shape[1]):
                    real = cube[0, fr, :, 0, :].real.astype(np.float64)
                    block = real.reshape(-1, BLOCK, real.shape[-1])
                    residual.append(np.sqrt(np.mean((block-reference)**2, axis=(1, 2))))
                    rms.append(np.sqrt(np.mean(real**2)))
                    repeat.append(np.mean(real.mean(axis=0)**2)/np.mean(real**2))
                    if fr % 25 == 0:
                        print(f'{shot}: read frame {fr}/{cube.shape[1]-1}', flush=True)
            residual = np.asarray(residual)
            threshold = np.percentile(residual[:marked], 99)
            normal_fraction = np.mean(residual <= threshold, axis=1)
            report = dict(shot=shot,frames=len(residual),marked_frame=marked,
                reference_rms_counts=float(np.sqrt(np.mean(reference**2))),
                baseline_bound_counts=float(threshold),
                baseline_residual_median_counts=float(np.median(residual[:marked])),
                last_frame_residual_median_counts=float(np.median(residual[-1])),
                last_frame_within_baseline_fraction=float(normal_fraction[-1]),
                post_frame_within_baseline_fraction=normal_fraction[marked+1:].tolist(),
                post_frames_all_blocks_within_baseline=(np.flatnonzero(normal_fraction[marked+1:]==1)+marked+1).tolist())
            g = product.create_group(name)
            g.attrs['source_file'] = p['file']
            g.attrs['source_id'] = previous.attrs['source_id']
            g.attrs['report_json'] = json.dumps(report)
            g.create_dataset('reference_real_counts', data=reference)
            g.create_dataset('block_residual_rms_counts', data=residual)
            g.create_dataset('frame_real_rms_counts', data=rms)
            g.create_dataset('frame_repeating_power_fraction', data=repeat)
            g.create_dataset('frame_fraction_within_baseline_bound', data=normal_fraction)
            x = np.arange(residual.size)/residual.shape[1] + .5/residual.shape[1]
            fig, ax = plt.subplots(figsize=(10, 4.4), layout='constrained')
            ax.plot(x, residual.ravel(), color='tab:blue', lw=.7)
            ax.axhline(threshold, color='tab:green', ls='--', label='Pre-event 99th percentile')
            ax.axvline(marked+onset/4096, color='black', ls='--', label='Marked signal change')
            ax.set(xlabel='Frame index (frame gaps omitted)', ylabel='RMS difference from pre-event waveform\n(ADC counts; 128-chirp blocks)',
                   title=f'Shot {shot} — RX0 real signal, entire recording', xlim=(0,len(residual)), ylim=(0,None))
            ax.legend(fontsize=9)
            fig.savefig(OUT / f'raw_recovery_{shot}.png', dpi=180)
            plt.close(fig)
            print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
