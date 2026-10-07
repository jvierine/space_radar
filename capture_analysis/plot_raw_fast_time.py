"""Plot original real receiver samples versus fast time and chirp index.

Run with conda run -n base python capture_analysis/plot_raw_fast_time.py.
No signal processing or voltage calibration is applied.
"""
from pathlib import Path
import argparse
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shot', choices=['6606', '6607'])
    frame_args = parser.add_mutually_exclusive_group()
    frame_args.add_argument('--frame', type=int, help='Zero-based frame index')
    frame_args.add_argument('--frames', type=int, nargs='+', help='Concatenate these zero-based frames in the supplied order')
    args = parser.parse_args()
    requested_frames = args.frames if args.frames is not None else ([args.frame] if args.frame is not None else None)
    if requested_frames is not None and args.shot is None:
        parser.error('--frame/--frames requires --shot')
    OUT.mkdir(exist_ok=True)
    selection_suffix = ('_frames_' + '_'.join(map(str, requested_frames)) if args.frames is not None
                        else f'_frame_{args.frame}' if args.frame is not None else '')
    suffix = f'_{args.shot}{selection_suffix}' if requested_frames is not None else ''
    with h5py.File(OUT / 'raw_receiver_iq.h5', 'r') as selection, \
         h5py.File(OUT / f'raw_fast_time{suffix}.h5', 'w') as product:
        product.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        product.attrs['processing'] = 'Real part of original RX0 samples; no subtraction, filtering, whitening, masking or conjugation'
        for name, shot in [('parallel', '6606'), ('perpendicular', '6607')]:
            if args.shot is not None and args.shot != shot:
                continue
            prior = selection[name]
            p = json.loads(prior.attrs['profile_json'])
            marked_frame = int(prior.attrs['frame'])
            frame_numbers = [marked_frame] if requested_frames is None else requested_frames
            real_blocks, hashes = [], {}
            onset = None
            with h5py.File(DATA / p['file'], 'r') as source:
                chirps_per_frame = source['radar_cube'].shape[2]
                for index, frame_number in enumerate(frame_numbers):
                    if not 0 <= frame_number < source['radar_cube'].shape[1]:
                        parser.error('Frame index outside the recording')
                    frame = source['radar_cube'][0, frame_number]
                    digest = hashlib.sha256(frame.tobytes()).hexdigest()
                    hashes[str(frame_number)] = digest
                    if frame_number == marked_frame:
                        assert digest == prior.attrs['source_frame_sha256']
                        onset = index * chirps_per_frame + int(prior.attrs['data_selected_onset_chirp'])
                    real_blocks.append(frame[:, 0, :].real.copy())
            real = np.concatenate(real_blocks, axis=0)
            n_chirps, n_samples = real.shape
            fast_us = np.arange(n_samples) / p['fs'] * 1e6
            # Edges centre each mesh cell on its actual sample/chirp coordinate.
            fast_edges = (np.arange(n_samples + 1) - .5) / p['fs'] * 1e6
            chirp_edges = np.arange(n_chirps + 1) - .5
            g = product.create_group(name)
            for attr in ['source_file', 'source_id']:
                g.attrs[attr] = prior.attrs[attr]
            g.attrs['source_frame_sha256_json'] = json.dumps(hashes)
            g.attrs['frames'] = frame_numbers
            g.attrs['frame_gap_note'] = 'Frames concatenated in supplied order; inter-frame gaps not represented'
            g.attrs['receiver'] = 0
            if onset is not None:
                g.attrs['data_selected_onset_chirp'] = onset
            g.attrs['fs_hz'] = p['fs']
            g.attrs['assumed_chirp_period_s'] = p['period']
            g.attrs['time_origin'] = f'First stored chirp of frame {frame_numbers[0]}'
            g.create_dataset('real_counts', data=real, compression='gzip').attrs['units'] = 'ADC counts; not calibrated volts'
            g.create_dataset('fast_time_us', data=fast_us).attrs['units'] = 'microseconds from first stored sample'
            g.create_dataset('chirp_index', data=np.arange(n_chirps))
            g.create_dataset('concatenated_time_s', data=np.arange(n_chirps) * p['period']).attrs['units'] = 'seconds; assumed chirp spacing, frame gaps omitted'
            g.create_dataset('source_frame_index', data=np.repeat(frame_numbers, chirps_per_frame))
            g.create_dataset('chirp_index_within_frame', data=np.tile(np.arange(chirps_per_frame), len(frame_numbers)))

            multiple = len(frame_numbers) > 1
            fig, ax = plt.subplots(figsize=(12, 7) if multiple else (8.8, 7.7), layout='constrained')
            mesh = ax.pcolormesh(chirp_edges, fast_edges, real.T, shading='flat',
                                 cmap='RdBu_r', vmin=-32768, vmax=32768, rasterized=True)
            if onset is not None:
                ax.axvline(onset, color='black', linestyle='--', linewidth=1)
            if multiple:
                for index, frame_number in enumerate(frame_numbers):
                    if index:
                        ax.axvline(index * chirps_per_frame - .5, color='black', linewidth=.9)
                    ax.text((index + .5) * chirps_per_frame, .98, f'Frame {frame_number}',
                            transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=10,
                            bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=2))
            title_frames = 'Frames ' + ', '.join(map(str, frame_numbers)) if multiple else f'Frame {frame_numbers[0]}'
            ax.set(xlabel='Concatenated chirp index (frame gaps omitted)' if multiple else 'Chirp index within frame', ylabel='Fast time within chirp (µs)',
                   xlim=(chirp_edges[0], chirp_edges[-1]), ylim=(fast_edges[0], fast_edges[-1]))
            fig.suptitle(f'Shot {shot} — RX0, real part of raw stored signal\n{title_frames}: all {n_chirps} chirps', fontsize=14)
            period = p['period']
            time_axis = ax.secondary_xaxis('top', functions=(lambda chirp: chirp * period,
                                                            lambda seconds: seconds / period))
            time_axis.set_xlabel(f'Time since start of frame {frame_numbers[0]} (s)')
            bar = fig.colorbar(mesh, ax=ax, pad=.03, fraction=.05)
            bar.set_label('Real signal amplitude (ADC counts)')
            bar.set_ticks([-32768, -16000, 0, 16000, 32768])
            fig.savefig(OUT / f'raw_fast_time_{shot}{selection_suffix}.png', dpi=200)
            plt.close(fig)
            print(f'{shot}: frames {frame_numbers}, RX0, {real.shape}, fast-time spacing {1e6/p["fs"]:.3f} us, onset chirp {onset}')


if __name__ == '__main__':
    main()
