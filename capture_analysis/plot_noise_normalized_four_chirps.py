"""Four-chirp Eq.25 energy / independent quiet-background output, in dB.

Run: conda run -n base python capture_analysis/plot_noise_normalized_four_chirps.py
The numerator is observed matched output energy, not signal-only energy.
The denominator includes residual clutter as well as receiver noise. No
calibrated thermal SNR, measured false-alarm threshold or pellet claim.
"""
from pathlib import Path
import hashlib
import json
import os
import time
import h5py
import numpy as np
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import search_quiet_frames as model

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / 'output'
SCRIPT = 'capture_analysis/plot_noise_normalized_four_chirps.py'
NC = 4
SEED = 660704
CONTEXT_START, CONTEXT_STOP = 1560, 2000


def db(x):
    return 10 * np.log10(np.maximum(x, 1e-12))


def main():
    begun = time.monotonic()
    product = OUT / 'chirp_spectra_6607_frames_34_35_36_37.h5'
    with h5py.File(product) as h:
        raw = h['complex_counts'][:4096].astype(np.complex128)
        source_id = str(h.attrs['source_id'])
        source_file = str(h.attrs['source_file'])
        # Stored scientific product remains usable when the source disk is absent.
        expected = json.loads(h.attrs['source_rx0_sha256_json'])['34']
        original = h['complex_counts'][:4096]
        # plot_chirp_spectra stored the source complex64 samples in complex128.
        np.testing.assert_array_equal(original, original.astype(np.complex64).astype(np.complex128))
        assert hashlib.sha256(original.astype(np.complex64).tobytes()).hexdigest() == expected
    with h5py.File(OUT / 'raw_receiver_iq.h5') as h:
        p = json.loads(h['perpendicular'].attrs['profile_json'])
    with h5py.File(OUT / 'background_subtracted_6607_frames_34_35_36_37.h5') as h:
        cutoff = int(h.attrs['quiet_chirp_stop_exclusive'])
    assert cutoff == 1624 and raw.shape == (4096, 128)
    # Disjoint baseline, training, validation and event samples. Do not fit the
    # baseline or noise reference using any of the four tested chirp trains.
    baseline = raw[:256].mean(axis=0)
    residual = raw - baseline
    train_starts = np.arange(256, 896, 8)       # 80 four-chirp blocks
    control_starts = np.arange(896, 1536, 8)   # 80 held-out four-chirp blocks
    event_starts = np.arange(cutoff - 16, cutoff + 16, NC)
    def blocks(starts):
        return np.stack([residual[s:s+NC].ravel() for s in starts]).T
    train, control, event = [blocks(s) for s in [train_starts, control_starts, event_starts]]
    assert train.shape == control.shape == (512, 80)
    assert control_starts[-1] + NC <= event_starts[0]
    model.NC = NC
    bank = object.__new__(model.QuietBank)
    bank.p = p
    u = p['T_adc'] + np.arange(128) / p['fs']
    bank.u = np.tile(u, NC)
    bank.t = (np.arange(NC)[:, None] * p['period'] + u).ravel()
    x = np.linspace(-1.2, 1.2, 49)
    vp = np.arange(0, 7000.1, 100)
    y = np.linspace(.05, 1, 20)
    theta = np.stack(np.meshgrid(x, vp, y, indexing='ij'), -1).reshape(-1, 3)
    # Proper-complex Gaussian surrogate with the uncentred second moment of
    # measured quiet residuals: C = train @ train.H / K. It is a model check,
    # distinct from actual held-out background, which need not be Gaussian.
    rng = np.random.default_rng(SEED)
    nsim = 128
    coeff = (rng.standard_normal((80, nsim)) + 1j*rng.standard_normal((80, nsim))) / np.sqrt(160)
    ratios = np.zeros((2, len(theta), len(event_starts)))
    noise_power = np.zeros((2, len(theta)))
    control_max = np.zeros(80)
    gaussian_max = np.zeros(nsim)
    fixed_index = int(np.argmin(np.sum((theta - [.6, 3500, .35])**2, axis=1)))
    fixed_gaussian = None
    eligible_count = 0
    with threadpool_limits(limits=8):
        for a in range(0, len(theta), 256):
            q = bank.template(theta[a:a+256])
            norms = np.sum(abs(q)**2, axis=1)
            valid = norms > .5
            eligible_count += int(valid.sum())
            np.testing.assert_allclose(norms[valid], 1, atol=1e-12)
            for sign in range(2):
                weights = q.conj() if sign == 0 else q
                training_output = weights @ train
                denominator = np.mean(abs(training_output)**2, axis=1)
                assert np.all(denominator[valid] > 0)
                def normalize(power):
                    return np.divide(power, denominator[:, None], out=np.zeros_like(power), where=valid[:, None])
                noise_power[sign, a:a+len(q)] = denominator
                ratios[sign, a:a+len(q)] = normalize(abs(weights @ event)**2)
                heldout = normalize(abs(weights @ control)**2)
                simulated = normalize(abs(training_output @ coeff)**2)
                control_max = np.maximum(control_max, heldout.max(axis=0))
                gaussian_max = np.maximum(gaussian_max, simulated.max(axis=0))
                # By construction, the mean training score for every valid
                # fixed filter is exactly one. This includes spectral colour
                # and correlations between samples/chirps in its denominator.
                np.testing.assert_allclose(normalize(abs(training_output)**2)[valid].mean(axis=1), 1, atol=1e-12)
                if sign == 0 and a <= fixed_index < a+len(q):
                    fixed_gaussian = simulated[fixed_index-a].copy()
            if a % 8192 == 0:
                print(f'Templates {a}/{len(theta)}; elapsed {time.monotonic()-begun:.1f}s', flush=True)
    canonical = ratios.reshape(2, len(x), len(vp), len(y), len(event_starts))
    signed = np.concatenate([canonical[:, ::-1, 1:, :, :][:, :, ::-1], canonical], axis=2)
    response = signed.max(axis=0)
    v = np.arange(-7000, 7000.1, 100)
    trial = np.array([[.3, 1400, .2], [-.9, 7000, .8]])
    np.testing.assert_allclose(bank.template(trial), bank.template(trial * [-1, -1, 1]), atol=1e-12, rtol=1e-12)
    assert np.isfinite(response).all() and fixed_gaussian is not None
    # Direct one-filter numerator/denominator check, independent of reshaping.
    ix0 = np.unravel_index(response[..., 0].argmax(), response[..., 0].shape)
    best = [x[ix0[0]], v[ix0[1]], y[ix0[2]]]
    q = bank.template(best)[0]
    direct = max(abs(np.vdot(w, event[:, 0]))**2 / np.mean(abs(w.conj() @ train)**2) for w in [q, q.conj()])
    np.testing.assert_allclose(response[ix0 + (0,)], direct, rtol=1e-12)
    # A separate white proper-complex noise check: the fixed-filter score is
    # exponential with linear mean 1; 0 dB refers to that mean, not mean dB.
    white = (rng.standard_normal(100000) + 1j*rng.standard_normal(100000)) / np.sqrt(2)
    assert abs(np.mean(abs(white)**2)-1) < .015
    assert abs(np.mean(abs(white)**2 > 3) - np.exp(-3)) < .003
    reports = []
    for i, start in enumerate(event_starts):
        ix = np.unravel_index(response[..., i].argmax(), response[..., i].shape)
        reports.append(dict(start_chirp=int(start), end_chirp=int(start+3),
                            start_time_s=float(start*p['period']),
                            end_sample_time_s=float((start+3)*p['period']+128/p['fs']),
                            peak_ratio=float(response[ix+(i,)]), peak_db=float(db(response[ix+(i,)])),
                            peak_theta=[float(x[ix[0]]), float(v[ix[1]]), float(y[ix[2]])],
                            disturbed=bool(start >= cutoff)))
    summary = dict(shot=6607, frame=34, cutoff_chirp=cutoff, baseline_chirps=[0,255],
                   training_chirp_interval=[256,895], heldout_interval=[896,1535],
                   training_blocks=80, heldout_blocks=80, gaussian_trials=nsim,
                   canonical_templates=len(theta), eligible_canonical_templates=eligible_count,
                   searched_iq_orientations=2, fixed_gaussian_linear_mean=float(fixed_gaussian.mean()),
                   gaussian_searched_peak_db_quantiles=np.percentile(db(gaussian_max), [50,90,95,100]).tolist(),
                   heldout_searched_peak_db_quantiles=np.percentile(db(control_max), [50,90,95,100]).tolist(),
                   event_reports=reports, elapsed_s=time.monotonic()-begun)
    output = OUT / 'noise_normalized_four_chirps_6607.h5'
    with h5py.File(output, 'w') as h:
        h.attrs['script'] = SCRIPT
        h.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        h.attrs['template_script_sha256'] = hashlib.sha256((ROOT/'search_quiet_frames.py').read_bytes()).hexdigest()
        h.attrs['input_product'] = str(product)
        h.attrs['input_sha256'] = hashlib.sha256(product.read_bytes()).hexdigest()
        h.attrs['source_frame34_rx0_sha256'] = expected
        h.attrs['source_file'] = source_file
        h.attrs['source_id'] = source_id
        h.attrs['profile_json'] = json.dumps(p)
        h.attrs['summary_json'] = json.dumps(summary)
        h.attrs['normalization'] = '|q.H z|^2 / mean_k |q.H n_k|^2, unit norm q; denominator separately estimated for every template and IQ orientation from disjoint quiet four-chirp blocks'
        h.attrs['interpretation'] = 'Observed matched output energy / quiet-background output energy. Includes thermal noise and structured clutter. Not signal-only SNR, calibrated thermal SNR or pellet detection. In stationary proper complex Gaussian noise, fixed-filter linear mean=1; searched maxima differ. Gaussian surrogate is a second-moment model, not proof of Gaussianity or stationarity.'
        h.attrs['baseline'] = 'Complex mean of frame34 chirps0-255 per fast-time bin, subtracted from training, held-out and event blocks. This independently estimated baseline differs from earlier full-quiet-interval plots.'
        h.attrs['scope'] = 'Eight fixed consecutive four-chirp trains1608-1639, four before and four after data-selected disturbance boundary1624; no time search for event selection. Max projections over omitted coordinate and stored IQ orientation. Cutoff1624 is not external shot time.'
        h.attrs['context_processing'] = 'Undecoded real part after the same complex quiet row-mean cancellation as the matched filter: mean of frame34chirps0-255 subtracted from every chirp. Context chirps1560-1999; seconds since frame34 start, no frame gap. Common symmetric count scale across all eight figures.'
        h.attrs['template_assumptions'] = 'Eq25 retarded phase,1/R^2 voltage,quasistatic175/350kHz HPFs,absIF<10MHz,>=80%valid; finite grid; unknown IQ sign. Noise-normalized phase matcher, not optimal covariance-whitened matched filter.'
        h.attrs['coordinate_definition'] = 'Projectile=(v0*t,0), radar=(x0,y0); time origin first chirp of EACH train. x0 is along-track coordinate, not slant range. Mirror(x0,v0)->(-x0,-v0) exactly ambiguous.'
        h.attrs['timing'] = p['timing_note']
        h.attrs['random_seed'] = SEED
        values = dict(x0_m=x, velocity_m_s=v, y0_m=y, window_starts=event_starts,
                      energy_to_quiet_noise_surface=response, canonical_quiet_output_power_counts2=noise_power,
                      max_ratio_v0_x0=response.max(axis=2), max_ratio_v0_y0=response.max(axis=0),
                      training_block_starts=train_starts, heldout_block_starts=control_starts,
                      heldout_max_ratio=control_max, gaussian_max_ratio=gaussian_max,
                      gaussian_fixed_ratio=fixed_gaussian, background_complex_counts=baseline,
                      event_residual_complex_counts=event.T.reshape(len(event_starts),4,128),
                      context_raw_complex_counts=raw[CONTEXT_START:CONTEXT_STOP],
                      context_residual_complex_counts=residual[CONTEXT_START:CONTEXT_STOP],
                      context_chirp_index=np.arange(CONTEXT_START,CONTEXT_STOP),
                      fast_time_us=np.arange(128)/p['fs']*1e6)
        for name, value in values.items():
            h.create_dataset(name, data=value, compression='gzip')
    print(json.dumps(summary, indent=2), flush=True)
    plot(summary, p, residual, response, x, v, y, control_max, gaussian_max)
    print(output, flush=True)


def provenance(fig):
    if os.getenv('SHOW_PROVENANCE', '1') != '0':
        fig.supxlabel('Script: '+SCRIPT, fontsize=8, color='.45')


def plot(summary, p, residual, response, x, v, y, controls, gaussian):
    plt.rcParams.update({'font.size':12, 'axes.titlesize':13})
    maxdb = max(db(response.max()), np.percentile(db(controls), 95), np.percentile(db(gaussian), 95))
    vmax = max(14, int(np.ceil(maxdb / 2)*2))
    context = residual[CONTEXT_START:CONTEXT_STOP].real.T
    count_limit = float(np.ceil(np.max(abs(context))/10000)*10000)
    for i, report in enumerate(summary['event_reports']):
        start = report['start_chirp']
        fig = plt.figure(figsize=(13.5,7.4), layout='constrained')
        gs = fig.add_gridspec(2,3, width_ratios=[1,1,.035], height_ratios=[.9,1.25])
        ax = fig.add_subplot(gs[0,:2])
        ix = np.arange(CONTEXT_START,CONTEXT_STOP+1)
        times = ix*p['period']
        fast_edges = np.arange(129)/p['fs']*1e6
        voltage=ax.pcolormesh(times,fast_edges,context, cmap='RdBu_r',vmin=-count_limit,vmax=count_limit,rasterized=True)
        voltage_bar=fig.add_subplot(gs[0,2])
        fig.colorbar(voltage,cax=voltage_bar,label='Real residual (ADC counts)')
        ax.axvspan(start*p['period'],(start+4)*p['period'],color='#ffd345',alpha=.25)
        for boundary in [start, start+4]:
            ax.axvline(boundary*p['period'],color='#996000',lw=2)
        ax.axvline(summary['cutoff_chirp']*p['period'],color='black',lw=1.3,ls='--')
        ax.text(.012,.93,'Real I/Q after complex quiet-mean cancellation; no decoding',transform=ax.transAxes,va='top',fontsize=10,bbox=dict(fc='white',ec='none',alpha=.9))
        ax.text(.985,.93,'Dashed: disturbed interval begins',transform=ax.transAxes,va='top',ha='right',fontsize=10,bbox=dict(fc='white',ec='none',alpha=.9))
        ax.set(xlabel='Time since start of frame 34 (s)',ylabel='Fast time (µs)',xlim=(times[0],times[-1]))
        ax.ticklabel_format(axis='x',style='plain',useOffset=False)
        top=ax.secondary_xaxis('top',functions=(lambda t:t/p['period'],lambda n:n*p['period']))
        top.set_xlabel('Chirp index within frame 34; gold lines select four chirps')
        a = fig.add_subplot(gs[1,0]); b = fig.add_subplot(gs[1,1])
        im=a.pcolormesh(v/1000,x,db(response[...,i].max(axis=2)),shading='nearest',cmap='magma',vmin=0,vmax=vmax,rasterized=True)
        b.pcolormesh(v/1000,y,db(response[...,i].max(axis=0)).T,shading='nearest',cmap='magma',vmin=0,vmax=vmax,rasterized=True)
        a.set(xlabel='Along-track velocity v₀ (km/s)',ylabel='Along-track position x₀ (m)',title='Maximum over perpendicular distance y₀')
        b.set(xlabel='Along-track velocity v₀ (km/s)',ylabel='Perpendicular distance y₀ (m)',title='Maximum over along-track position x₀')
        ca=fig.add_subplot(gs[1,2]);fig.colorbar(im,cax=ca,label='Matched energy / noise (dB)')
        fig.suptitle(f'Shot 6607 — frame 34, chirps {start}–{start+3}; peak {report["peak_db"]:.1f} dB',fontsize=17)
        provenance(fig)
        name=OUT/f'noise_match_6607_chirps_{start}_{start+3}'
        fig.savefig(str(name)+'.png',dpi=200);fig.savefig(str(name)+'.pdf');plt.close(fig)
    fig, ax = plt.subplots(figsize=(11.5,5.5),layout='constrained')
    for values,label,color in [(gaussian,'Gaussian surrogate (128 trains)', '#306fbe'),(controls,'Held-out quiet data (80 trains)','#d56a25')]:
        ordered=np.sort(db(values))
        ax.step(ordered,np.arange(1,len(ordered)+1)/len(ordered),where='post',label=label,color=color,lw=2)
    quiet_reports=[r for r in summary['event_reports'] if not r['disturbed']]
    for i, report in enumerate(quiet_reports):
        ax.axvline(report['peak_db'],color='#5e9653',ls=['-', '--', ':', '-.'][i],alpha=.8,label=f'Chirps {report["start_chirp"]}–{report["end_chirp"]}: {report["peak_db"]:.1f} dB')
    ax.axvline(13,color='.35',lw=1,ls=':',label='13 dB reference; not a detection threshold')
    ax.set(xlabel='Largest noise-normalized response over the full bank (dB)',ylabel='Fraction of control trains below this value',ylim=(0,1.02),title='The search also produces peaks in quiet data')
    ax.legend(fontsize=11,loc='lower right');ax.grid(alpha=.2)
    provenance(fig);fig.savefig(OUT/'four_chirp_noise_controls_6607.png',dpi=200);fig.savefig(OUT/'four_chirp_noise_controls_6607.pdf');plt.close(fig)


if __name__ == '__main__':
    main()
