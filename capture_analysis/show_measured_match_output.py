"""Display the measured-data matched-filter output and conditional range/IF.

Run with base conda Python. This sweeps the previously fitted, fixed eight-
chirp template over every chirp lag in the marked frame. It does not repeat
the full geometry bank at every lag, and is not a detection threshold.
"""
from pathlib import Path
import hashlib
import json

import h5py
import numpy as np
from scipy.signal import correlate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from analyze_captures import DATA, echo
from apply_measured_filter import prepare
from search_simulated_echoes import Bank

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent/'output'


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.4), layout='constrained')
    reports = {}
    with h5py.File(ROOT/'memo_004_measured_filter.h5') as previous, \
         h5py.File(OUT/'measured_match_output.h5', 'w') as out:
        out.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        out.attrs['input_sha256'] = hashlib.sha256((ROOT/'memo_004_measured_filter.h5').read_bytes()).hexdigest()
        out.attrs['scope'] = 'Fixed previously fitted template, all chirp lags; no new full-bank search or confirmed detection'
        for row, name in enumerate(['parallel', 'perpendicular']):
            group = previous[name]
            p = json.loads(group.attrs['profile_json'])
            roi = int(group.attrs['frame_of_interest'])
            saved = group[f'frame_{roi}']
            report = json.loads(saved.attrs['report_json'])
            start = report['start_chirp']
            theta = report['recovered_theta']
            mean, good, psd = [group[key][:] for key in
                ['clutter_mean_counts', 'valid_adc_mask', 'background_psd_counts2']]
            with h5py.File(DATA/p['file']) as source:
                raw = source['radar_cube'][0, roi]
            assert hashlib.sha256(raw.tobytes()).hexdigest() == saved.attrs['frame_sha256']
            white = prepare(raw, mean, good, psd)
            h = Bank(p, np.zeros((8, p['samples']))).template(theta).reshape(8, p['samples'])
            sign = -1 if report['iq_conjugated'] else 1
            if sign < 0:
                h = h.conj()
            score = np.zeros(len(raw)-7)
            direct = 0.
            for rx in range(4):
                q = np.fft.fft(h*good[rx], axis=-1, norm='ortho')/np.sqrt(psd[rx])
                norm = np.vdot(q, q).real
                inner = correlate(white[:, rx], q, mode='valid', method='fft')[:, 0]
                score += abs(inner)**2/norm
                direct += abs(np.vdot(q, white[start:start+8, rx]))**2/norm
            np.testing.assert_allclose(score[start], direct, rtol=1e-10)
            np.testing.assert_allclose(direct, report['score'], rtol=1e-10)
            event = echo(p, theta, 8)
            retained = abs(event['iq']) > 0
            fstored = sign*event['if_hz']
            rmin, rmax = event['range'][retained].min(), event['range'][retained].max()
            fmin, fmax = fstored[retained].min(), fstored[retained].max()
            residual = saved['clutter_subtracted_counts'][:]*good[None]
            spectrum = np.mean(abs(np.fft.fft(residual, axis=-1, norm='ortho'))**2, axis=(0, 1))
            frequencies = np.fft.fftfreq(p['samples'], 1/p['fs'])
            peakfreq = frequencies[np.argmax(spectrum)]
            shot = '6606' if name == 'parallel' else '6607'
            a = axes[row, 0]
            times = np.arange(len(score))*p['period']*1e3
            a.plot(times, score, lw=.8, label='Fixed-template lag sweep')
            a.scatter(saved['window_starts'][:]*p['period']*1e3,
                      saved['window_coarse_max'][:], s=22, facecolors='none',
                      edgecolors='tab:orange', label='Earlier bank maxima (tested windows)')
            a.scatter([times[start]], [score[start]], color='tab:red', marker='*', s=100, label='Window fitted in Memo 4')
            a.axvline(report['onset_chirp']*p['period']*1e3, color='k', ls='--', lw=.8, label='Data-change onset')
            a.set(yscale='log', xlabel='Window start within marked frame (ms)',
                  ylabel='Match statistic Λ', title=f'Shot {shot}: measured match output')
            if row == 0:
                a.legend(fontsize=8, loc='upper left')
            a = axes[row, 1]
            t = event['t']*1e6
            for chirp in range(8):
                a.plot(t[chirp], event['range'][chirp], color='.65', lw=1)
                a.plot(t[chirp], np.where(retained[chirp], event['range'][chirp], np.nan), color='tab:orange', lw=2)
            a.plot([], [], color='.65', label='Conditional fitted trajectory')
            a.plot([], [], color='tab:orange', lw=2, label='Samples retained by template')
            a.set(xlabel='Time within fitted eight-chirp window (µs)', ylabel='Model range at reflection (m)',
                  title=f'Shot {shot}: fitted range; speed at bound')
            a.legend(fontsize=8)
            a = axes[row, 2]
            order = np.argsort(frequencies)
            a.plot(frequencies[order]/1e6, 10*np.log10(np.maximum(spectrum[order]/spectrum.max(), 1e-12)), color='k', lw=1, label='Measured, clutter-subtracted')
            a.axvspan(fmin/1e6, fmax/1e6, color='tab:orange', alpha=.25, label='Retained fitted-template IF span')
            a.axvline(peakfreq/1e6, color='tab:blue', ls=':', label=f'Largest FFT bin: {peakfreq/1e6:.3f} MHz')
            a.set(xlabel='Frequency in stored I/Q (MHz)', ylabel='Measured spectral power / peak (dB)',
                  title=f'Shot {shot}: same fitted window', xlim=(-6.25, 6.25), ylim=(-45, 2))
            a.legend(fontsize=8, loc='lower left')
            result = dict(shot=shot, frame=roi, window_start_chirp=start, score=direct,
                range_retained_m=[float(rmin), float(rmax)],
                fitted_stored_if_retained_MHz=[float(fmin/1e6), float(fmax/1e6)],
                largest_measured_fft_bin_MHz=float(peakfreq/1e6),
                fitted_speed_m_s=theta[1], retained_samples_per_rx=int(retained.sum()),
                fixed_template_sweep_max=float(score.max()),
                fixed_template_sweep_max_chirp=int(score.argmax()))
            reports[name] = result
            g = out.create_group(name)
            g.attrs['report_json'] = json.dumps(result)
            g.attrs['source_frame_sha256'] = saved.attrs['frame_sha256']
            for key, values, units in [('lag_score', score, '1'), ('window_start_time', times/1e3, 's'),
                ('template_sample_time', event['t'], 's'), ('template_range', event['range'], 'm'),
                ('template_stored_if', fstored, 'Hz'), ('template_retained_mask', retained, '1'),
                ('measured_frequency', frequencies, 'Hz'), ('measured_spectral_power', spectrum, 'ADC counts squared')]:
                ds = g.create_dataset(key, data=values, compression='gzip'); ds.attrs['units'] = units
            print(name, json.dumps(result), flush=True)
        out.attrs['reports_json'] = json.dumps(reports)
    fig.suptitle('Measured chamber data: matched-filter response and conditional fitted range / IF\n'
                 'Left: one fixed fitted template swept over time; right panels: Memo 4 fitted window. No confirmed projectile echo.', fontsize=12)
    fig.savefig(OUT/'measured_match_output.png', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
