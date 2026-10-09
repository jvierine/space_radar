"""Ideal local information for Memo 7; unknown common complex amplitude.
Run: conda run -n base python capture_analysis/memo_007_information.py
"""
from pathlib import Path
import h5py
import numpy as np

C = 299792458.0
F0, S, FS, PERIOD = 77e9, 99.987e12, 12.5e6, 25.37e-6
THETA = np.array([0.1, 341.51, 0.15])  # x0 [m], v0 [m/s], y0 [m]


def phase(theta, t, u):
    x, v, y = theta
    r = np.hypot(x-v*t, y)
    tau = 2*r/C
    return 2*np.pi*(-F0*tau-S*u*tau+0.5*S*tau*tau)


def information(n):
    u0 = 2e-6 + np.arange(225)/FS
    t = (np.arange(n)[:, None]*PERIOD+u0).ravel()
    u = np.tile(u0, n)
    x, v, y = THETA
    X = x-v*t
    R = np.hypot(X, y)
    jac = (2*np.pi*(-F0-S*u+2*S*R/C)*2/C)[:, None]*np.array([X/R, -t*X/R, y/R]).T
    # Numerical derivative regression; phase differences, not wrapped phase.
    for k, eps in enumerate([1e-6, 1e-3, 1e-6]):
        delta = np.zeros(3)
        delta[k] = eps
        fd = (phase(THETA+delta, t, u)-phase(THETA-delta, t, u))/(2*eps)
        np.testing.assert_allclose(jac[:, k], fd, rtol=3e-6, atol=2e-6)
    jac -= jac.mean(axis=0)  # eliminate unknown complex phase/amplitude
    rho = 100*n  # 20 dB integrated energy SNR for one chirp
    _, singular, vt = np.linalg.svd(jac, full_matrices=False)
    covariance = (vt.T/singular**2)@vt * len(t)/(2*rho)
    sd = np.sqrt(np.diag(covariance))
    correlation = covariance/np.outer(sd, sd)
    return t, u, jac, covariance, sd, correlation


if __name__ == '__main__':
    output = Path(__file__).with_suffix('.h5')
    with h5py.File(output, 'w') as h:
        h.attrs['model'] = 'Receive-time range; phase coherent across chirps; unit amplitude; white circular complex Gaussian noise'
        h.attrs['parameter_order'] = 'x0 [m], v0 [m/s], y0 [m]'
        h.attrs['rho_single_chirp'] = 100.0
        h.attrs['chirp_period_assumed_s'] = PERIOD
        h.create_dataset('truth', data=THETA)
        for n in [1, 2, 4, 8, 16]:
            t, u, jac, cov, sd, corr = information(n)
            g = h.create_group(str(n))
            for key, value in [('t', t), ('u', u), ('phase_jacobian_centered', jac), ('covariance', cov), ('standard_deviation', sd), ('correlation', corr)]:
                g.create_dataset(key, data=value)
            print(n, 'std(x0,v0,y0):', sd, 'corr(xv,xy,vy):', corr[np.triu_indices(3, 1)])
    print('Derivative checks passed; saved', output)
