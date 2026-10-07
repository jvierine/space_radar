"""Simple measured four-chirp Eq.25 response for the FMCW lecture.

Run in base conda. Fixed windows immediately before the data-selected
disturbance; no selection of a stronger window elsewhere and no detection
claim. Reuse the exact retarded-phase template from search_quiet_frames.
"""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
from threadpoolctl import threadpool_limits
import search_quiet_frames as model

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / 'output'
NC = 4


def main():
    product = OUT / 'background_subtracted_6607_frames_34_35_36_37.h5'
    with h5py.File(product) as h:
        cutoff = int(h.attrs['quiet_chirp_stop_exclusive'])
        assert cutoff == 1624
        starts = np.arange(cutoff - 4 * NC, cutoff, NC)
        data = np.stack([h['residual_complex_counts'][s:s+NC] for s in starts])
        background = h['background_complex_counts'][:]
        source_id = h.attrs['source_id']
    with h5py.File(OUT / 'raw_receiver_iq.h5') as h:
        p = json.loads(h['perpendicular'].attrs['profile_json'])
    assert data.shape == (4, 4, 128)
    energy = np.sum(abs(data)**2, axis=(1, 2))
    # This is a separate process; change the reused template's chirp count
    # explicitly rather than pretending an eight-chirp result is four chirps.
    model.NC = NC
    bank = object.__new__(model.QuietBank)
    bank.p = p
    u = p['T_adc'] + np.arange(128) / p['fs']
    bank.u = np.tile(u, NC)
    bank.t = (np.arange(NC)[:, None] * p['period'] + u).ravel()
    x = np.linspace(-1.2, 1.2, 49)
    vpos = np.arange(0, 7000.1, 100)
    y = np.linspace(.05, 1, 20)
    theta = np.stack(np.meshgrid(x, vpos, y, indexing='ij'), -1).reshape(-1, 3)
    powers = np.empty((2, len(theta), len(starts)))
    d = data.reshape(4, -1).T
    with threadpool_limits(limits=8):
        for a in range(0, len(theta), 256):
            q = bank.template(theta[a:a+256])
            assert q.shape[1] == 512
            norms = np.sum(abs(q)**2, axis=1)
            np.testing.assert_allclose(norms[norms > 0], 1, atol=1e-12)
            powers[0, a:a+len(q)] = abs(q.conj() @ d)**2 / energy
            powers[1, a:a+len(q)] = abs(q @ d)**2 / energy
    canonical = powers.reshape(2, len(x), len(vpos), len(y), len(starts))
    # Monostatic range history is invariant under (x0,v)->(-x0,-v).
    surface = np.concatenate([canonical[:, ::-1, 1:, :, :][:, :, ::-1], canonical], axis=2)
    v = np.arange(-7000, 7000.1, 100)
    response = surface.max(axis=0)
    assert np.isfinite(response).all() and 0 <= response.min() and response.max() <= 1 + 1e-12
    # Check the mirror shortcut against directly generated negative templates.
    trial = np.array([[.3, 1400, .2], [-.9, 7000, .8]])
    np.testing.assert_allclose(bank.template(trial), bank.template(trial * [-1, -1, 1]), atol=1e-12, rtol=1e-12)
    last = response[..., -1]
    ix = np.unravel_index(last.argmax(), last.shape)
    best = np.array([x[ix[0]], v[ix[1]], y[ix[2]]])
    sign = int(surface[:, ix[0], ix[1], ix[2], -1].argmax())
    q = bank.template(best)[0]
    if sign:
        q = q.conj()
    fitted = q * np.vdot(q, d[:, -1])
    np.testing.assert_allclose(np.sum(abs(fitted)**2) / energy[-1], last[ix], atol=1e-12)
    # An exact bank waveform has unit coherence; check the score definition.
    injected = bank.template(np.array([.6, 3500, .35]))[0]
    np.testing.assert_allclose(abs(np.vdot(injected, injected))**2, 1, atol=1e-12)
    report = dict(shot=6607, frame=34, chirps=[int(starts[-1]), cutoff-1],
                  cutoff_chirp=cutoff, theta=best.tolist(), mirror_theta=(best*[-1, -1, 1]).tolist(),
                  matched_fraction=float(last[ix]), iq_conjugated=bool(sign),
                  window_span_us=float(((NC-1)*p['period']+128/p['fs'])*1e6),
                  integration_time_us=float(NC*128/p['fs']*1e6))
    output = OUT / 'lecture_four_chirp_match_6607.h5'
    with h5py.File(output, 'w') as h:
        h.attrs['script'] = 'capture_analysis/lecture_four_chirp_match.py'
        h.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        h.attrs['template_source'] = 'capture_analysis/search_quiet_frames.py: QuietBank.template'
        h.attrs['template_source_sha256'] = hashlib.sha256((ROOT/'search_quiet_frames.py').read_bytes()).hexdigest()
        h.attrs['input_product'] = str(product)
        h.attrs['input_sha256'] = hashlib.sha256(product.read_bytes()).hexdigest()
        h.attrs['source_id'] = source_id
        h.attrs['profile_json'] = json.dumps(p)
        h.attrs['report_json'] = json.dumps(report)
        h.attrs['normalization'] = '|sum(z*conj(q))|^2 / (sum|z|^2 * sum|q|^2); fraction of residual energy explained by one template with fitted complex amplitude'
        h.attrs['scope'] = 'Four fixed non-overlapping four-chirp windows 1608–1623 in frame34; main slide uses last window1620–1623. Complex quiet mean0–1623 removed from all samples, as in displayed spectra. All windows stop before data-selected cutoff1624; no independent trigger.'
        h.attrs['assumptions'] = 'Same Eq25 retarded phase,1/R^2 amplitude,175/350kHz HPFs,absIF<10MHz,>=80%valid samples as previous measured bank. Both IQ orientations. Finite x0[-1.2,1.2]m,y0[.05,1]m,v[-7,7]km/s. Unknown absolute geometry and IQ orientation; no whitening or calibrated SNR.'
        h.attrs['interpretation'] = 'Template resemblance of measured quiet residual, not pellet detection, recovered projectile parameters, thermal SNR or false-alarm probability. Signed velocity mirror is exactly ambiguous.'
        for key, value in [('x0_m', x), ('velocity_m_s', v), ('y0_m', y),
                           ('window_starts', starts), ('match_fraction_surface', response),
                           ('profile_x0', last.max(axis=(1, 2))), ('profile_velocity', last.max(axis=(0, 2))),
                           ('profile_y0', last.max(axis=(0, 1))), ('residual_complex_counts', data),
                           ('fitted_complex_counts_last_window', fitted.reshape(4, 128)),
                           ('background_complex_counts', background), ('time_s', bank.t.reshape(4,128))]:
            h.create_dataset(key, data=value, compression='gzip')
    print(json.dumps(report, indent=2), flush=True)
    print(output, flush=True)


if __name__ == '__main__':
    main()
