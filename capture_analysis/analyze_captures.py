"""Recorded-profile Eq.25 echo and profiled matched-filter feasibility study.

Run: conda run -n base python capture_analysis/analyze_captures.py
No measured echo is claimed recovered. Geometry/gains/material are assumptions.
"""
from pathlib import Path
import hashlib
import os
import json
import h5py
import numpy as np
from scipy.special import spherical_jn, spherical_yn
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C = 299792458.0
KB = 1.380649e-23
ROOT = Path(__file__).resolve().parent
DATA = Path(os.environ.get('FMCW_CAPTURE_DIR', ROOT.parent / 'data/raw'))
PT = 10**((12-30)/10)  # datasheet maximum, NOT measured TX setting
G = 10**(6/10)  # assumed 6 dBi each TX and RX
TSYS = 290*10**(15/10)  # antenna at 290 K plus 15 dB receiver NF
COUNTS = np.array([1, 2, 4, 8, 16, 32])


def sphere(frequency, diameter, extra=0):
    """Complex PEC-sphere backscatter length; |length|^2 = monostatic RCS.

    Riccati-Bessel partial waves. Global scattering phase is a nuisance.
    """
    k = 2*np.pi*np.asarray(frequency)/C
    x = k*diameter/2
    n = np.arange(1, int(np.ceil(np.max(x)+4*np.max(x)**(1/3)+12))+extra+1)
    xx = np.asarray(x)[..., None]
    j = spherical_jn(n, xx)
    y = spherical_yn(n, xx)
    psi = xx*j
    xi = xx*(j+1j*y)
    dp = j+xx*spherical_jn(n, xx, derivative=True)
    dx = j+1j*y+xx*(spherical_jn(n, xx, derivative=True)+1j*spherical_yn(n, xx, derivative=True))
    a, b = -dp/dx, -psi/xi
    return np.sqrt(np.pi)/k*np.sum((2*n+1)*(-1.)**n*(a-b), axis=-1)


def read_profile(path):
    with h5py.File(path) as h:
        p = {name: float(h[name][()]) for name in
             ['fs', 'freq_slope', 'f_start', 'T_adc', 'T_idle', 'if_max', 'speed', 'projectile_size']}
        p['samples'] = h['radar_cube'].shape[-1]
        p['comment'] = str(h.attrs.get('comment', ''))
    p['file'] = path.name
    p['period'] = p['samples']/p['fs']+p['T_idle']+p['T_adc']+0.37e-6
    p['timing_note'] = '0.37 us ramp excess inferred from SpaceRadar reader; verify rampEndTime'
    return p


def echo(p, theta, chirps=32, gate=True, phase_only=False):
    """Eq.25: global target time, ramp-local modulation time; exact light time.

    theta=(x0,v0,y0), R(t)=sqrt((v0*t-x0)^2+y0^2).
    Beam/illumination is approximated by projectile x in [0,1] m.
    Receiver: quasistatic two first-order HPFs and ideal Complex-1x side mask.
    This is NOT a measured transient/DFE response or noise covariance model.
    """
    x0, v, y0 = theta
    u = p['T_adc']+np.arange(p['samples'])/p['fs']
    t = np.arange(chirps)[:, None]*p['period']+u[None, :]
    u = np.broadcast_to(u, t.shape)
    x = v*t-x0
    rr = x*x+y0*y0
    d = rr/(np.sqrt((x*v)**2+(C*C-v*v)*rr)+x*v)
    tau = 2*d
    xs = x-v*d
    r = np.hypot(xs, y0)
    vr = v*xs/r
    td = 2*vr/(C+vr)
    s, f0 = p['freq_slope'], p['f_start']
    phase = 2*np.pi*(-f0*tau-s*u*tau+0.5*s*tau*tau)
    f = -s*tau-(f0+s*(u-tau))*td
    hp = (1j*f/(175e3+1j*f))*(1j*f/(350e3+1j*f))
    side = (f < 0) & (f > -p['if_max']) & (u >= tau)
    active = ((v*(t-d) >= 0) & (v*(t-d) <= 1)) if gate else np.ones_like(t, dtype=bool)
    if phase_only:
        z = np.exp(1j*phase)*hp*side*active
    else:
        rf = f0+s*(u-tau)
        scattering = sphere(rf, p['projectile_size']*1e-3)
        z = np.sqrt(PT*G*G)*(C/rf)/(4*np.pi)**1.5*scattering/r**2*np.exp(1j*phase)*hp*side*active
    return dict(t=t, u=u, range=r, radial_velocity=vr, if_hz=f, iq=z,
                phase=phase, active=active, retained=side)


def score(data, template):
    e = np.vdot(template, template).real
    return 0.0 if e == 0 else abs(np.vdot(template, data))**2/e


def local_metric(p, theta, chirps=8):
    # Normalized complex derivative, projected off unknown complex amplitude.
    h = echo(p, theta, chirps, gate=False, phase_only=True)['iq'].ravel()
    h /= np.linalg.norm(h)
    derivatives = []
    for i, eps in enumerate([1e-6, 0.01, 1e-6]):
        a, b = np.array(theta, float), np.array(theta, float)
        a[i] += eps; b[i] -= eps
        ha = echo(p, a, chirps, False, True)['iq'].ravel()
        hb = echo(p, b, chirps, False, True)['iq'].ravel()
        dh = (ha/np.linalg.norm(ha)-hb/np.linalg.norm(hb))/(2*eps)
        derivatives.append(dh-h*np.vdot(h, dh))
    d = np.array(derivatives)
    metric = np.real(d.conj()@d.T)
    # Grid edges allocate 10% total local power mismatch over 3 axes.
    steps = 2*np.sqrt(0.1/(3*np.diag(metric)))
    scaled = metric*steps[:, None]*steps[None, :]
    return metric, steps, np.linalg.eigvalsh(scaled)


def main():
    ROOT.mkdir(exist_ok=True)
    paths = sorted(DATA.glob('2025*.nc'))
    assert len(paths) == 2
    profiles = [read_profile(path) for path in paths]
    # Validate PEC calculation against Rayleigh limit and truncation refinement.
    a, f = 1e-5, 1e9
    exact = abs(sphere(f, 2*a))**2
    rayleigh = 9*np.pi*a*a*(2*np.pi*f*a/C)**4
    assert abs(exact/rayleigh-1) < 1e-5
    for p in profiles:
        d = p['projectile_size']*1e-3
        assert np.allclose(sphere(p['f_start'], d), sphere(p['f_start'], d, 10), atol=1e-12)
    summary = {'assumptions': {'tx_power_dbm': 12, 'tx_gain_dbi': 6,
        'rx_gain_dbi': 6, 'system_temperature_K': TSYS,
        'sphere': 'PEC; metadata size assumed diameter', 'extra_loss_db': 0,
        'noise': 'white complex input-equivalent variance k*Tsys*fs; no clutter',
        'receiver': 'assumed Complex1x, 175/350kHz HPF; quasistatic transfer'}, 'profiles': profiles}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    with h5py.File(ROOT/'feasibility.h5', 'w') as h:
        h.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        h.attrs['assumptions_json'] = json.dumps(summary['assumptions'])
        for i, p in enumerate(profiles):
            name = 'parallel' if i == 0 else 'perpendicular'
            # Illustrative geometries, NOT inferred surveyed coordinates.
            theta = [-1.0132, p['speed'], 0.1138] if i == 0 else [0.4732, p['speed'], 0.2138]
            e = echo(p, theta)
            wavelength = C/p['f_start']
            sigma = abs(sphere(p['f_start'], p['projectile_size']*1e-3))**2
            pr1m = PT*G*G*wavelength**2*sigma/(4*np.pi)**3
            ideal_snr = pr1m*COUNTS*p['samples']/p['fs']/(KB*TSYS)
            energy = np.sum(abs(e['iq'])**2, axis=1)
            # Best contiguous window for each N, not necessarily starting at launch.
            moving_snr = np.array([np.max(np.convolve(energy, np.ones(n), 'full'))/p['fs']/(KB*TSYS) for n in COUNTS])
            p.update(sphere_rcs_m2=float(sigma), reference_power_1m_W=float(pr1m),
                     reference_snr_db=(10*np.log10(ideal_snr)).tolist(),
                     moving_snr_db=(10*np.log10(moving_snr)).tolist(), assumed_theta=theta)
            gr = h.create_group(name)
            gr.attrs['profile_json'] = json.dumps(p)
            for key, value in e.items(): gr.create_dataset(key, data=value)
            gr.create_dataset('chirp_counts', data=COUNTS)
            gr.create_dataset('ideal_constant_1m_snr', data=ideal_snr)
            gr.create_dataset('moving_1m_path_snr', data=moving_snr)
            t = e['t'].ravel()*1e6
            axes[i, 0].scatter(t, e['if_hz'].ravel()/1e6, c=e['active'].ravel(), s=2, cmap='coolwarm')
            axes[i, 0].axhspan(-10, 0, alpha=.1, color='green')
            axes[i, 0].set(xlabel='Time (µs)', ylabel='Signed IF (MHz)', title=name+'; Eq.25 rx × LO*')
            axes[i, 1].plot(COUNTS, 10*np.log10(ideal_snr), 'o-', label='Constant 1m range, ideal')
            axes[i, 1].plot(COUNTS, 10*np.log10(moving_snr), 's-', label='Moving, 1m path, assumed geometry')
            axes[i, 1].set(xscale='log', xlabel='Integrated chirps', ylabel='Thermal MF SNR (dB)', title='One RX; assumed 6 dBi TX/RX')
            axes[i, 1].legend(fontsize=8)
        p = profiles[1]
        theta = np.array([0.5, p['speed'], 0.2])
        metric, steps, eigen = local_metric(p, theta)
        rng = np.random.default_rng(6607)
        truth = theta+steps*np.array([.23, -.17, .31])
        model = echo(p, truth, 8, gate=False, phase_only=True)['iq'].ravel()
        model *= np.sqrt(100/np.vdot(model, model).real)  # total injected SNR 20 dB
        noise = (rng.normal(size=model.size)+1j*rng.normal(size=model.size))/np.sqrt(2)
        data = model+noise
        grids = [theta[i]+steps[i]*np.arange(-4, 5) for i in range(3)]
        likelihood = np.empty((9, 9, 9))
        for ix, x in enumerate(grids[0]):
            for iv, v in enumerate(grids[1]):
                for iy, y in enumerate(grids[2]):
                    q = echo(p, (x, v, y), 8, False, True)['iq'].ravel()
                    likelihood[ix, iv, iy] = score(data, q)
        idx = np.unravel_index(np.argmax(likelihood), likelihood.shape)
        recovered = [float(grids[i][idx[i]]) for i in range(3)]
        summary['injection_test'] = dict(truth=truth.tolist(), recovered=recovered,
            grid_steps=steps.tolist(), metric=metric.tolist(), scaled_metric_eigenvalues=eigen.tolist(),
            injected_snr_db=20, detected_statistic=float(likelihood[idx]), noise_only_mean=1,
            scope='Local 729-template recovery, not global identifiability or measured detection')
        recovered_template = echo(p, recovered, 8, False, True)['iq'].ravel()
        coherence = abs(np.vdot(recovered_template, model))**2/(np.vdot(recovered_template, recovered_template).real*np.vdot(model, model).real)
        summary['injection_test']['recovered_template_power_coherence'] = float(coherence)
        assert coherence > 0.9, 'Synthetic echo should be detected with little mismatch'
        assert abs(echo(p, theta, 8, False, True)['iq']).max() > 0
        gr = h.create_group('injection_test')
        for name, val in [('data', data), ('template_signal', model), ('metric', metric),
                          ('grid_steps', steps), ('scores', likelihood), ('truth', truth), ('recovered', recovered)]:
            gr.create_dataset(name, data=val)
        for name, val in zip(['x0_grid', 'v0_grid', 'y0_grid'], grids): gr.create_dataset(name, data=val)
    fig.suptitle('Recorded-profile simulation: assumed geometry/material/gains, not calibrated detection\nScript: capture_analysis/analyze_captures.py', fontsize=11)
    fig.savefig(ROOT/'feasibility.png', dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    ax[0].plot(np.arange(128)/p['fs']*1e6, data.reshape(8, 128)[-1].real, label='Noisy I')
    ax[0].plot(np.arange(128)/p['fs']*1e6, model.reshape(8, 128)[-1].real, label='Injected I')
    ax[0].set(xlabel='ADC time within last chirp (µs)', ylabel='Normalized sample amplitude', title='Synthetic 20 dB total-SNR injection')
    ax[0].legend(fontsize=8)
    im = ax[1].imshow(np.max(likelihood, axis=2), origin='lower', aspect='auto',
        extent=[grids[1][0], grids[1][-1], grids[0][0], grids[0][-1]])
    ax[1].set(xlabel='v0 (m/s)', ylabel='x0 (m)', title='Matched-filter score, maximized over y0')
    fig.colorbar(im, ax=ax[1], label='Score; noise-only mean = 1')
    fig.suptitle('Local recovery test; ridge demonstrates parameter coupling\nScript: capture_analysis/analyze_captures.py', fontsize=10)
    fig.savefig(ROOT/'matched_filter.png', dpi=160)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__': main()
