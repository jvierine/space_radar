"""Two maximum-projection match-score maps for the four-chirp lecture.

Run: conda run -n base python capture_analysis/plot_four_chirp_marginals.py
SHOW_PROVENANCE=0 hides generating-script text; default is visible.
Each map maximizes over the omitted coordinate, as requested by the user.
"""
from pathlib import Path
import hashlib
import json
import os
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / 'output'
SCRIPT = 'capture_analysis/plot_four_chirp_marginals.py'


def main():
    source = OUT / 'lecture_four_chirp_match_6607.h5'
    with h5py.File(source) as h:
        score = h['match_fraction_surface'][:, :, :, -1]
        x = h['x0_m'][:]
        v = h['velocity_m_s'][:] / 1000
        y = h['y0_m'][:]
        report = json.loads(h.attrs['report_json'])
        assert report['chirps'] == [1620, 1623] and score.shape == (49, 141, 20)
    # Maximum over the third spatial coordinate, preserving score units.
    # I/Q orientation has already been maximized in the input product.
    vx = score.max(axis=2)
    vy = score.max(axis=0).T
    assert vx.shape == (len(x), len(v)) and vy.shape == (len(y), len(v))
    assert np.isfinite(vx).all() and np.isfinite(vy).all()
    np.testing.assert_allclose(vx.max(), score.max(), rtol=1e-12)
    np.testing.assert_allclose(vy.max(), score.max(), rtol=1e-12)
    assert np.all(vx[:, :, None] >= score) and np.all(vy.T[None, :, :] >= score)
    np.testing.assert_allclose(vx, vx[::-1, ::-1], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(vy, vy[:, ::-1], rtol=1e-12, atol=1e-12)
    product = OUT / 'four_chirp_marginals_6607.h5'
    with h5py.File(product, 'w') as h:
        h.attrs['script'] = SCRIPT
        h.attrs['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        h.attrs['input_product'] = str(source)
        h.attrs['input_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
        h.attrs['report_json'] = json.dumps(report)
        h.attrs['projection'] = 'Maximum normalized squared coherence over y0 for (v0,x0), over x0 for (v0,y0); IQ orientation remains maximized as in input. Profiled match score, not likelihood/posterior.'
        h.attrs['omitted_coordinate_bounds'] = 'y0 .05..1m in .05m increments; x0 -1.2..1.2m in .05m increments'
        h.attrs['normalization'] = 'Match fractions 0..1, displayed as percentages; same linear colour limits for both surfaces'
        h.attrs['interpretation'] = 'Shot6607 frame34 chirps1620–1623, same quiet complex mean. Maximum residual-template resemblance along omitted coordinate; no confirmed pellet echo, calibrated SNR, or posterior probability.'
        for name, data in [('x0_m', x), ('velocity_m_s', v*1000), ('y0_m', y),
                           ('max_match_v0_x0', vx), ('max_match_v0_y0', vy)]:
            h.create_dataset(name, data=data)
    bg, fg, muted = '#07111f', '#f4f7fa', '#9db0c3'
    plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm',
                         'font.size': 14, 'text.color': fg, 'axes.labelcolor': fg,
                         'xtick.color': muted, 'ytick.color': muted, 'axes.edgecolor': muted,
                         'axes.facecolor': bg, 'figure.facecolor': bg, 'savefig.facecolor': bg})
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout='constrained')
    ceiling = float(np.ceil(100 * max(vx.max(), vy.max())))
    for ax, data, vertical, label, title in zip(axes, [vx, vy], [x, y],
            [r'Along-path position $x_0$ (m)', r'Perpendicular distance $y_0$ (m)'],
            [r'$v_0,x_0$: maximum over $y_0$', r'$v_0,y_0$: maximum over $x_0$']):
        mesh = ax.pcolormesh(v, vertical, 100*data, shading='nearest',
                             cmap='magma', vmin=0, vmax=ceiling, rasterized=True)
        ax.set(xlabel=r'Along-path velocity $v_0$ (km/s)', ylabel=label, title=title,
               xlim=(-7, 7), ylim=(vertical[0], vertical[-1]), xticks=[-7, -3.5, 0, 3.5, 7])
        ax.tick_params(top=False, right=False)
    bar = fig.colorbar(mesh, ax=axes, pad=.025, shrink=.88)
    bar.set_label('Maximum matched residual energy (%)')
    if os.environ.get('SHOW_PROVENANCE', '1') != '0':
        fig.supxlabel('Script: '+SCRIPT, fontsize=9, color=muted)
    for ext in ['png', 'pdf']:
        fig.savefig(OUT/f'four_chirp_marginals_6607.{ext}', dpi=200)
    plt.close(fig)
    print(json.dumps(dict(product=str(product), v0_x0_peak_percent=float(100*vx.max()),
                          v0_y0_peak_percent=float(100*vy.max()), colour_max_percent=ceiling)), flush=True)


if __name__ == '__main__':
    main()
