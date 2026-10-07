"""Real residual and one complex spectrum per chirp after quiet-background removal.

Run: conda run -n base python capture_analysis/plot_background_subtracted_views.py
Uses the unmodified complex samples in the existing per-chirp spectrum products.
Preserves their frequency sign, time coordinates and spectral colour reference.
"""
from pathlib import Path
import hashlib
import json

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / 'output'


def decorate(ax, frames, chirps_per_frame, onset, period, spectral=False):
    colour = 'white' if spectral else 'black'
    for i, frame in enumerate(frames):
        if i:
            ax.axvline(i * chirps_per_frame - .5, color=colour, lw=.8, alpha=.8)
        ax.text((i + .5) * chirps_per_frame, .98, f'Frame {frame}',
                transform=ax.get_xaxis_transform(), ha='center', va='top',
                fontsize=10, bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=2))
    # Same unlabeled data-selected boundary as the raw views; not a trigger time.
    ax.axvline(onset, color=colour, ls='--', lw=.9, alpha=.8)
    ax.set_xlabel('Concatenated chirp index (frame gaps omitted)')
    ax.set_xlim(-.5, len(frames) * chirps_per_frame - .5)
    top = ax.secondary_xaxis('top', functions=(lambda x: x * period, lambda t: t / period))
    top.set_xlabel(f'Time since start of frame {frames[0]} (s; frame gaps omitted)')


def main():
    script = 'capture_analysis/plot_background_subtracted_views.py'
    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for shot, frames, name in [('6606', [43, 44, 45, 46], 'parallel'),
                               ('6607', [34, 35, 36, 37], 'perpendicular')]:
        suffix = f'{shot}_frames_' + '_'.join(map(str, frames))
        source_product = OUT / f'chirp_spectra_{suffix}.h5'
        with h5py.File(source_product) as source:
            z = source['complex_counts'][:]
            frequency = source['frequency_hz'][:]
            period = float(source.attrs['chirp_period_s_assumed'])
            fs = float(source.attrs['fs_hz'])
            reference = float(source.attrs['reference_spectral_power_counts2'])
            provenance = dict(source.attrs)
            np.testing.assert_array_equal(source['source_frame_index'][:], np.repeat(frames, 4096))
        with h5py.File(OUT / 'raw_receiver_iq.h5') as meta:
            onset = int(meta[name].attrs['data_selected_onset_chirp'])
            np.testing.assert_equal(meta[name].attrs['frame'], frames[0])
        # The 6607 audit identified a sustained rise five chirps before the
        # original marker. Exclude those from this new background estimate.
        quiet_stop = onset
        cutoff_source = 'raw_receiver_iq.h5: data_selected_onset_chirp'
        if shot == '6607':
            with h5py.File(OUT / 'quiet_cutoff_audit_6607.h5') as audit:
                quiet_stop = int(audit.attrs['disturbance_first_chirp'])
            cutoff_source = 'quiet_cutoff_audit_6607.h5: disturbance_first_chirp'
        assert z.shape == (16384, 128) and 0 < quiet_stop <= onset < 4096
        background = z[:quiet_stop].mean(axis=0)
        residual = z - background[None, :]
        window = np.hanning(z.shape[1])
        transform = np.fft.fftshift(np.fft.fft(residual * window, axis=1), axes=1) / window.sum()
        power = abs(transform)**2
        db = 10 * np.log10(np.maximum(power, np.finfo(float).tiny) / reference)

        # Verify fixed complex subtraction over both quiet and disturbed data,
        # FFT linearity (not power subtraction), and window-normalized Parseval.
        np.testing.assert_allclose(residual[:quiet_stop].mean(axis=0), 0, atol=1e-10)
        np.testing.assert_allclose(residual + background, z, rtol=0, atol=1e-10)
        raw_fft = np.fft.fftshift(np.fft.fft(z * window, axis=1), axes=1) / window.sum()
        background_fft = np.fft.fftshift(np.fft.fft(background * window)) / window.sum()
        np.testing.assert_allclose(transform, raw_fft - background_fft, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(power.sum(axis=1), len(frequency) * np.sum(abs(residual * window)**2, axis=1) / window.sum()**2, rtol=1e-12)
        assert np.isfinite(db).all()

        product = OUT / f'background_subtracted_{suffix}.h5'
        with h5py.File(product, 'w') as out:
            for key, value in provenance.items():
                out.attrs[key] = value
            out.attrs['input_product_script_sha256'] = provenance['script_sha256']
            out.attrs['input_product'] = str(source_product)
            out.attrs['script'] = script
            out.attrs['script_sha256'] = script_hash
            out.attrs['quiet_frame'] = frames[0]
            out.attrs['quiet_chirp_start_inclusive'] = 0
            out.attrs['quiet_chirp_stop_exclusive'] = quiet_stop
            out.attrs['quiet_cutoff_source'] = cutoff_source
            out.attrs['unlabeled_change_marker_chirp'] = onset
            out.attrs['processing'] = 'Fixed complex mean per fast-time bin over quiet chirps only, subtracted from ALL chirps; real residual display and separate Hann FFT128 per complex chirp; no whitening, conjugation or range conversion'
            out.attrs['spectral_reference'] = 'Maximum raw spectral power in the corresponding unsubtracted four-frame plot; same -70 to 0 dB colour limits'
            out.attrs['real_display_limit_counts'] = 32768
            out.attrs['quiet_residual_complex_rms_counts'] = float(np.sqrt(np.mean(abs(residual[:quiet_stop])**2)))
            out.attrs['verification'] = 'Quiet residual mean zero; reconstruction; complex FFT linearity; window-normalized Parseval checked for every chirp'
            out.create_dataset('background_complex_counts', data=background)
            out.create_dataset('residual_complex_counts', data=residual, compression='gzip')
            out.create_dataset('spectral_power_counts2', data=power, compression='gzip')
            out.create_dataset('spectral_power_db_relative_to_raw_peak', data=db.astype(np.float32), compression='gzip')
            out.create_dataset('frequency_hz', data=frequency)
            out.create_dataset('fast_time_us', data=np.arange(128) / fs * 1e6)
            out.create_dataset('chirp_index', data=np.arange(len(z)))
            out.create_dataset('concatenated_time_s', data=np.arange(len(z)) * period)
            out.create_dataset('source_frame_index', data=np.repeat(frames, 4096))
        subtitle = f'Frames {", ".join(map(str, frames))}: all {len(z)} chirps; quiet mean from frame {frames[0]}, chirps 0–{quiet_stop - 1}'
        edges = np.arange(len(z) + 1) - .5
        fig, ax = plt.subplots(figsize=(12, 7))
        fig.subplots_adjust(left=.085, right=.90, bottom=.13, top=.80)
        fast_edges = (np.arange(129) - .5) / fs * 1e6
        mesh = ax.pcolormesh(edges, fast_edges, residual.real.T, shading='flat', cmap='RdBu_r', vmin=-32768, vmax=32768, rasterized=True)
        decorate(ax, frames, 4096, onset, period)
        ax.set_ylabel('Fast time within chirp (µs)')
        ax.set_ylim(fast_edges[0], fast_edges[-1])
        fig.suptitle(f'Shot {shot} — RX0, real component after complex-background subtraction\n' + subtitle, fontsize=13)
        bar = fig.colorbar(mesh, ax=ax, pad=.03, fraction=.05)
        bar.set_label('Real residual amplitude (ADC counts)')
        bar.set_ticks([-32768, -16000, 0, 16000, 32768])
        fig.text(.5, .012, f'Script: {script} · Same quiet background removed from every chirp; no decoding', ha='center', fontsize=7, color='gray')
        real_image = OUT / f'background_subtracted_real_{suffix}.png'
        fig.savefig(real_image, dpi=200)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(12, 7))
        fig.subplots_adjust(left=.085, right=.90, bottom=.13, top=.80)
        frequency_edges = (np.arange(129) - 64 - .5) * fs / 128 / 1e6
        mesh = ax.pcolormesh(edges, frequency_edges, db.T, shading='flat', cmap='magma', vmin=-70, vmax=0, rasterized=True)
        decorate(ax, frames, 4096, onset, period, spectral=True)
        ax.set_ylabel('Signed sampled beat frequency (MHz)')
        ax.set_ylim(-fs / 2e6, fs / 2e6)
        fig.suptitle(f'Shot {shot} — RX0, spectrum of each background-subtracted chirp\n' + subtitle, fontsize=13)
        bar = fig.colorbar(mesh, ax=ax, pad=.03, fraction=.05)
        bar.set_label('Spectral power (dB relative to raw-plot maximum)')
        fig.text(.5, .012, f'Script: {script} · Complex residual; Hann FFT per chirp; same spectral colour reference as raw plot', ha='center', fontsize=7, color='gray')
        spectrum_image = OUT / f'background_subtracted_spectra_{suffix}.png'
        fig.savefig(spectrum_image, dpi=200)
        plt.close(fig)
        print(json.dumps(dict(shot=shot, quiet_stop_exclusive=quiet_stop, quiet_duration_s=quiet_stop * period,
                              quiet_residual_rms_counts=float(np.sqrt(np.mean(abs(residual[:quiet_stop])**2))),
                              real_image=str(real_image), spectrum_image=str(spectrum_image), product=str(product))), flush=True)


if __name__ == '__main__':
    main()
