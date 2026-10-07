"""Bounded coarse-to-fine 3-D search of noisy, radar-equation scaled echoes.

Run with base conda Python. Synthetic receiver/noise, not measured background.
Coarse bank then bounded local refinement of separated bank maxima.
"""
from pathlib import Path
import hashlib
import json
import time
import h5py
import numpy as np
from scipy.optimize import minimize
from analyze_captures import C, KB, TSYS, PT, G, echo, sphere

ROOT = Path(__file__).resolve().parent


class Bank:
    def __init__(self, p, data, first_chirp=0):
        self.p = p
        self.data = np.asarray(data).ravel()
        nc = len(self.data)//p['samples']
        u = p['T_adc']+np.arange(p['samples'])/p['fs']
        self.t = ((first_chirp+np.arange(nc))[:, None]*p['period']+u[None, :]).ravel()
        self.u = np.tile(u, nc)
        rf = p['f_start']+p['freq_slope']*self.u
        # Delay perturbs the sweep frequency by <= a few hundred kHz in this
        # metre-scale example. Cache receive-time scattering; verify against
        # the exact delayed-frequency forward model at injected parameters.
        self.amplitude_rf = np.sqrt(PT*G*G)*(C/rf)/(4*np.pi)**1.5*sphere(rf, p['projectile_size']*1e-3)

    def template(self, theta):
        theta = np.atleast_2d(theta)
        x0, v, y = [theta[:, i, None] for i in range(3)]
        x = v*self.t-x0
        rr = x*x+y*y
        d = rr/(np.sqrt((x*v)**2+(C*C-v*v)*rr)+x*v)
        tau = 2*d
        xs = x-v*d
        r = np.hypot(xs, y)
        vr = v*xs/r
        td = 2*vr/(C+vr)
        s, f0 = self.p['freq_slope'], self.p['f_start']
        f = -s*tau-(f0+s*(self.u-tau))*td
        phase = 2*np.pi*(-f0*tau-s*self.u*tau+.5*s*tau*tau)
        hp = (1j*f/(175e3+1j*f))*(1j*f/(350e3+1j*f))
        active = (v*(self.t-d) >= 0) & (v*(self.t-d) <= 1)
        retained = (f < 0) & (f > -self.p['if_max']) & (self.u >= tau)
        return self.amplitude_rf/r**2*np.exp(1j*phase)*hp*active*retained

    def scores(self, theta):
        h = self.template(theta)
        energy = np.sum(abs(h)**2, axis=1)
        inner = np.sum(h.conj()*self.data, axis=1)
        return np.divide(abs(inner)**2, energy, out=np.zeros_like(energy), where=energy > 0)


def search(p, data, grids, first_chirp=0):
    bank = Bank(p, data, first_chirp)
    parameters = np.stack(np.meshgrid(*grids, indexing='ij'), axis=-1).reshape(-1, 3)
    values = np.empty(len(parameters))
    for start in range(0, len(parameters), 64):
        values[start:start+64] = bank.scores(parameters[start:start+64])
    bounds = np.array([[g[0], g[-1]] for g in grids])
    origin = bounds.mean(axis=1)
    scale = np.array([.005, 20, .005])
    unit_bounds = (bounds-origin[:, None])/scale[:, None]
    coarse_step = np.array([g[1]-g[0] for g in grids])
    seeds = []
    for idx in np.argsort(values)[::-1]:
        candidate = parameters[idx]
        if all(np.linalg.norm((candidate-other)/coarse_step) >= 3 for other in seeds):
            seeds.append(candidate)
        if len(seeds) == 12: break
    results = []
    for seed in seeds:
        fit = minimize(lambda w: -bank.scores(origin+scale*w)[0],
            (seed-origin)/scale, method='Nelder-Mead', bounds=unit_bounds,
            options={'maxiter': 600, 'xatol': 1e-6, 'fatol': 1e-6})
        results.append((float(-fit.fun), origin+scale*fit.x, bool(fit.success)))
    results.sort(key=lambda a: a[0], reverse=True)
    return bank, values.reshape(*(len(g) for g in grids)), results


def main():
    start = time.monotonic()
    rng = np.random.default_rng(20261006)
    reports = {}
    with h5py.File(ROOT/'feasibility.h5') as inp, h5py.File(ROOT/'search_experiment.h5', 'w') as out:
        out.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        out.attrs['status'] = 'Synthetic experiment; thermal white noise and assumed receiver/gains; no real IQ used'
        out.attrs['template_approximation'] = 'receive-time RF scattering cached; exact delayed phase/range retained'
        for name in ['parallel', 'perpendicular']:
            p = json.loads(inp[name].attrs['profile_json'])
            if name == 'parallel':
                truth = np.array([-1.0132, p['speed'], .1138])
                grids = [np.linspace(-1.10, -.90, 41), np.linspace(6400, 7000, 31), np.linspace(.05, .20, 31)]
            else:
                truth = np.array([.4732, p['speed'], .2138])
                grids = [np.linspace(.40, .55, 31), np.linspace(6000, 6600, 31), np.linspace(.15, .30, 31)]
            exact = echo(p, truth, 32)['iq']
            noise_var = KB*TSYS*p['fs']
            noise = np.sqrt(noise_var/2)*(rng.normal(size=exact.shape)+1j*rng.normal(size=exact.shape))
            data = exact+noise
            gr = out.create_group(name)
            gr.attrs['profile_json'] = json.dumps(p)
            for key, value in [('truth', truth), ('exact_echo', exact),
                    ('noisy_iq', data), ('noise_variance_W', noise_var)]:
                gr.create_dataset(key, data=value, compression='gzip' if np.ndim(value) else None)
            for key, value in zip(['x0_grid', 'v0_grid', 'y0_grid'], grids): gr.create_dataset(key, data=value)
            # A fixed known window: launch for receding; chirp 6 for passing.
            # The latter avoids pretending an approaching, rejected sideband
            # contains an observable echo. It is NOT a blind time search.
            first_chirp = 0 if name == 'parallel' else 6
            reports[name] = {}
            for n in [1, 2, 4, 8]:
                observed = data[first_chirp:first_chirp+n]
                hexact = exact[first_chirp:first_chirp+n].ravel()
                bank, coarse, results = search(p, (observed/np.sqrt(noise_var)).ravel(), grids, first_chirp)
                best_score, recovered, success = results[0]
                htrue = bank.template(truth).ravel()
                cache_error = np.linalg.norm(htrue-hexact)/np.linalg.norm(hexact)
                assert cache_error < .002
                wave = bank.template(recovered).ravel()
                coherence = abs(np.vdot(wave, hexact))**2/(np.vdot(wave, wave).real*np.vdot(hexact, hexact).real)
                fitted_scale = np.vdot(wave, observed.ravel())/np.vdot(wave, wave)
                thermal_snr = np.sum(abs(hexact)**2)/noise_var
                w = hexact/np.sqrt(noise_var)
                truth_score = abs(np.vdot(w, observed.ravel()/np.sqrt(noise_var)))**2/np.vdot(w, w).real
                report = dict(truth=truth.tolist(), recovered=recovered.tolist(),
                    chirps=n, first_chirp=first_chirp, score=best_score, expected_noise_score=1,
                    success=success, power_coherence=float(coherence), coarse_templates=int(coarse.size),
                    thermal_snr_db=float(10*np.log10(thermal_snr)), truth_template_score=float(truth_score),
                    cached_scattering_relative_error=float(cache_error),
                    grid_bounds=[[float(g[0]), float(g[-1])] for g in grids],
                    finite_search_note='Coarse bank + 12 local refinements; not guaranteed global maximum or confidence interval')
                reports[name][n] = report
                case = gr.create_group(f'case_{n}')
                case.attrs['report_json'] = json.dumps(report)
                for key, value in [('recovered', recovered), ('fitted_echo', (fitted_scale*wave).reshape(n,128)),
                        ('coarse_scores', coarse), ('thermal_snr', thermal_snr),
                        ('refined_parameters', np.array([r[1] for r in results])),
                        ('refined_scores', np.array([r[0] for r in results]))]:
                    case.create_dataset(key, data=value, compression='gzip' if np.ndim(value) else None)
                print(name, n, json.dumps(report), flush=True)
        out.attrs['elapsed_seconds'] = time.monotonic()-start
    print('Elapsed seconds:', time.monotonic()-start)


if __name__ == '__main__': main()
