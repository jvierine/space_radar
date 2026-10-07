"""Exploratory Eq.25 matched filter on measured chamber I/Q, not synthetic data.

Uses independent pre-shot training, identical energy-window selection and
bounded template searches on event/control/post-shot frames. Frame gaps are
never joined. Conditional result: ADC mode, geometry and ramp timing unverified.
"""
from pathlib import Path
import hashlib
import json
import time

import h5py
import numpy as np
from scipy.optimize import minimize

from analyze_captures import DATA, read_profile
from search_simulated_echoes import Bank

ROOT = Path(__file__).resolve().parent
NC = 8


def rail_mask(z):
    return (z.real <= -32768) | (z.real >= 32767) | (z.imag <= -32768) | (z.imag >= 32767)


def prepare(z, mean, good, psd):
    residual = (z-mean[None])*good[None]
    return np.fft.fft(residual, axis=-1, norm='ortho')/np.sqrt(psd)[None]


def windows(white, fs):
    # The same selection rule is applied to all event/control/post-shot frames.
    high = abs(np.fft.fftfreq(white.shape[-1], 1/fs)) >= 1e6
    energy = np.sum(abs(white[..., high])**2, axis=(1, 2))
    smooth = np.convolve(energy, np.ones(NC)/NC, mode='valid')
    starts = []
    for block in range(0, len(energy), 512):
        end = min(block+512, len(smooth))
        if end > block:
            starts.append(block+int(np.argmax(smooth[block:end])))
    baseline = max(float(np.median(energy[:512])), 1e-12)
    high_idx = np.flatnonzero(smooth > 10*baseline)
    onset = int(high_idx[0]) if len(high_idx) else -1
    if onset >= 0:
        starts.extend([max(0, onset-NC), onset, min(len(smooth)-1, onset+NC)])
    return np.array(sorted(set(starts))), energy, onset


class MeasuredBank:
    def __init__(self, p, good, psd):
        self.p, self.good, self.psd = p, good, psd
        self.base = Bank(p, np.zeros((NC, p['samples']), dtype=complex))
        self.grids = [np.linspace(-1.2, 1.2, 49),
                      np.linspace(.95*p['speed'], 1.05*p['speed'], 17),
                      np.linspace(.05, .5, 19)]
        self.parameters = np.stack(np.meshgrid(*self.grids, indexing='ij'), -1).reshape(-1, 3)
        self.bounds = np.array([[g[0], g[-1]] for g in self.grids])
        self.step = np.array([g[1]-g[0] for g in self.grids])
        b, length = len(self.parameters), NC*p['samples']
        self.normalized = np.empty((2, good.shape[0], b, length), dtype=np.complex64)
        for first in range(0, b, 256):
            theta = self.parameters[first:first+256]
            h = self.base.template(theta).reshape(-1, NC, p['samples'])
            for sign in range(2):
                wave = h if sign == 0 else h.conj()
                for rx in range(good.shape[0]):
                    f = np.fft.fft(wave*good[rx], axis=-1, norm='ortho')/np.sqrt(psd[rx])
                    q = f.reshape(len(theta), -1)
                    norm = np.sqrt(np.sum(abs(q)**2, axis=1))
                    q = np.divide(q, norm[:, None], out=np.zeros_like(q), where=norm[:, None] > 0)
                    self.normalized[sign, rx, first:first+len(theta)] = q

    def coarse(self, white, starts):
        d = np.stack([white[i:i+NC] for i in starts]).transpose(2, 0, 1, 3)
        # d shape [RX, candidate window, chirp, ADC sample].
        out = np.zeros((2, len(self.parameters), len(starts)))
        for sign in range(2):
            for rx in range(len(self.good)):
                observed = d[rx].reshape(len(starts), -1).T.astype(np.complex64)
                inner = self.normalized[sign, rx].conj() @ observed
                out[sign] += abs(inner)**2
        return out

    def evaluate(self, theta, white, sign, full=False):
        h = self.base.template(theta).reshape(NC, self.p['samples'])
        if sign: h = h.conj()
        score, alphas = 0., []
        for rx in range(len(self.good)):
            q = np.fft.fft(h*self.good[rx], axis=-1, norm='ortho')/np.sqrt(self.psd[rx])
            en = np.vdot(q, q).real
            inner = np.vdot(q, white[:, rx])
            score += abs(inner)**2/en if en > 0 else 0
            alphas.append(inner/en if en > 0 else 0j)
        fitted = np.asarray(alphas)[None, :, None]*h[:, None, :]
        return (float(score), fitted, np.asarray(alphas)) if full else float(score)

    def refine(self, white, coarse, window_index):
        origin = self.bounds.mean(axis=1)
        scale = np.array([.025, 35., .025])
        ub = (self.bounds-origin[:, None])/scale[:, None]
        candidates = []
        for sign in range(2):
            values = coarse[sign, :, window_index]
            seeds = []
            for idx in np.argsort(values)[::-1]:
                theta = self.parameters[idx]
                if all(np.linalg.norm((theta-other)/self.step) >= 3 for other in seeds):
                    seeds.append(theta)
                if len(seeds) == 8: break
            for seed in seeds:
                fit = minimize(lambda w: -self.evaluate(origin+scale*w, white, sign),
                    (seed-origin)/scale, method='Nelder-Mead', bounds=ub,
                    options={'maxiter': 350, 'xatol': 1e-5, 'fatol': 1e-5})
                candidates.append((float(-fit.fun), origin+scale*fit.x, sign, bool(fit.success)))
        return sorted(candidates, key=lambda r:r[0], reverse=True)


def main():
    began = time.monotonic()
    outpath = ROOT/'memo_004_measured_filter.h5'
    with h5py.File(outpath, 'w') as out:
        out.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        out.attrs['status'] = 'Actual measured counts; bounded exploratory search, not confirmed projectile detection'
        out.attrs['assumptions'] = 'Complex1x retained side; minimum HPFs; reader-derived 13.71us period; no surveyed geometry'
        out.attrs['method'] = 'Pre-shot mean and diagonal FFT noise whitening; separate complex gain per RX; test both I/Q signs'
        for name, file in zip(['parallel', 'perpendicular'], sorted(DATA.glob('2025*.nc'))):
            p = read_profile(file)
            with h5py.File(file, 'r') as raw:
                roi = int(raw['frame_of_interest'][()])
                ds = raw['radar_cube']
                train_idx = roi-4
                train = ds[0, train_idx].astype(np.complex128)
                mean = train.mean(axis=0)
                rail_rate = rail_mask(train).mean(axis=0)
                good = rail_rate < .01
                good[:, 0] = False  # model excludes previous-ramp return at ADC start
                tf = np.fft.fft((train-mean[None])*good[None], axis=-1, norm='ortho')
                psd = np.mean(abs(tf)**2, axis=0)
                psd = np.maximum(psd, .01*np.median(psd, axis=-1, keepdims=True))
                group = out.create_group(name)
                group.attrs['profile_json'] = json.dumps(p)
                group.attrs['source_file'] = file.name
                group.attrs['source_id'] = str(raw.attrs['id'])
                group.attrs['source_bytes'] = file.stat().st_size
                group.attrs['source_mtime_ns'] = file.stat().st_mtime_ns
                group.attrs['source_shape'] = ds.shape
                group.attrs['frame_of_interest'] = roi
                group.attrs['training_frame'] = train_idx
                group.attrs['training_frame_sha256'] = hashlib.sha256(train.astype(np.complex64).tobytes()).hexdigest()
                for key, value in [('clutter_mean_counts',mean),('background_psd_counts2',psd),
                                   ('valid_adc_mask',good),('rail_rate',rail_rate)]:
                    group.create_dataset(key,data=value)
                print(name, 'building', 49*17*19, 'templates; valid ADC cells',int(good.sum()),'/',good.size,flush=True)
                bank = MeasuredBank(p, good, psd)
                for key, grid in zip(['x0_grid_m','v0_grid_m_s','y0_grid_m'],bank.grids): group.create_dataset(key,data=grid)
                reports = []
                for fr in [roi-3,roi-2,roi-1,roi,roi+1]:
                    original = ds[0, fr]
                    z = original.astype(np.complex128)
                    white = prepare(z, mean, good, psd)
                    starts, energy, onset = windows(white,p['fs'])
                    coarse = bank.coarse(white, starts)
                    best_index = np.unravel_index(np.argmax(coarse),coarse.shape)
                    wi = best_index[2]
                    start = int(starts[wi])
                    selected = white[start:start+NC]
                    results = bank.refine(selected,coarse,wi)
                    best, theta, sign, success = results[0]
                    _, fitted, alphas = bank.evaluate(theta,selected,sign,full=True)
                    observed_energy = float(np.sum(abs(selected)**2))
                    role = 'event' if fr==roi else 'postshot' if fr>roi else 'control'
                    report = dict(frame=fr,role=role,start_chirp=start,chirps=NC,
                        score=best,coarse_score=float(coarse[best_index]),
                        template_explained_whitened_fraction=best/observed_energy,
                        recovered_theta=theta.tolist(),iq_conjugated=bool(sign),optimizer_success=success,
                        selected_windows=len(starts),onset_chirp=onset,
                        raw_rail_fraction=float(rail_mask(z).mean()),
                        median_band_energy=float(np.median(energy)),
                        grid_bounds=bank.bounds.tolist(),
                        search_scope='15827 coarse templates, two I/Q signs; refine 8 seeds/sign in the top coarse window; finite search')
                    report['near_search_boundary'] = bool(np.any(np.min(abs(theta[:, None]-bank.bounds),axis=1)/bank.step < .05))
                    reports.append(report)
                    gr=group.create_group(f'frame_{fr}')
                    gr.attrs['report_json']=json.dumps(report)
                    gr.attrs['frame_sha256']=hashlib.sha256(original.tobytes()).hexdigest()
                    for k,v in [('band_energy',energy),('window_starts',starts),
                                ('window_coarse_max',coarse.max(axis=(0,1))),
                                ('observed_counts',z[start:start+NC]),
                                ('clutter_subtracted_counts',z[start:start+NC]-mean[None]),
                                ('fitted_counts',fitted),('complex_gains',alphas),
                                ('refined_scores',np.array([r[0] for r in results])),
                                ('refined_theta',np.array([r[1] for r in results]))]:
                        gr.create_dataset(k,data=v,compression='gzip')
                    # Save the actual event-frame spectrum for a measured-data diagnostic.
                    if role=='event':
                        gr.create_dataset('whitened_spectrogram',data=np.mean(abs(white)**2,axis=1),compression='gzip')
                        if onset >= 0:
                            oi = int(np.flatnonzero(starts==onset)[0])
                            onset_data = white[onset:onset+NC]
                            onset_results = bank.refine(onset_data,coarse,oi)
                            oscore, otheta, osign, osuccess = onset_results[0]
                            _, ofit, _ = bank.evaluate(otheta,onset_data,osign,full=True)
                            orep = dict(frame=fr,start_chirp=onset,chirps=NC,score=oscore,
                                recovered_theta=otheta.tolist(),iq_conjugated=bool(osign),optimizer_success=osuccess,
                                template_explained_whitened_fraction=oscore/float(np.sum(abs(onset_data)**2)),
                                near_search_boundary=bool(np.any(np.min(abs(otheta[:, None]-bank.bounds),axis=1)/bank.step < .05)),
                                scope='Window at data-selected change onset; not independently known launch time')
                            og=gr.create_group('onset_fit')
                            og.attrs['report_json']=json.dumps(orep)
                            og.create_dataset('observed_counts',data=z[onset:onset+NC],compression='gzip')
                            og.create_dataset('clutter_subtracted_counts',data=z[onset:onset+NC]-mean[None],compression='gzip')
                            og.create_dataset('fitted_counts',data=ofit,compression='gzip')
                            print(name,'ONSET',json.dumps(orep),flush=True)
                    print(name,json.dumps(report),flush=True)
                    del original,z,white,coarse
                group.attrs['reports_json']=json.dumps(reports)
                event=next(r for r in reports if r['role']=='event')
                control=max(r['score'] for r in reports if r['role']=='control')
                group.attrs['event_to_max_control_score']=event['score']/control
                group.attrs['interpretation']='Event/control contrast is conditional, not a calibrated false-alarm probability. Persistent post-shot change and clipping require diagnosis.'
                del bank
        out.attrs['elapsed_seconds']=time.monotonic()-began
        def units(name, ds):
            if isinstance(ds,h5py.Dataset):
                ds.attrs['units'] = ('ADC counts' if 'counts' in name and 'psd' not in name else
                    'ADC counts squared' if 'psd' in name else 'm/s' if name.endswith('v0_grid_m_s') else
                    'm' if name.endswith('_grid_m') else 'chirp index' if name.endswith('window_starts') else
                    'columns: m, m/s, m' if name.endswith('refined_theta') else '1')
        out.visititems(units)
    print('Wrote',outpath,'seconds',time.monotonic()-began,flush=True)


if __name__ == '__main__': main()
