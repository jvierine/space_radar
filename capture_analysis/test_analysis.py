"""Engineering/model checks, not hardware or measured-echo validation."""
import unittest
from pathlib import Path
import h5py
import numpy as np
from analyze_captures import C, echo, score, sphere


class ModelTests(unittest.TestCase):
    def setUp(self):
        self.p = dict(fs=12.5e6, freq_slope=60e12, f_start=79.5e9,
                      T_adc=0, if_max=10e6, period=13.71e-6,
                      samples=128, projectile_size=6)

    def test_rayleigh_and_truncation(self):
        f, a = 1e9, 1e-5
        expected = 9*np.pi*a*a*(2*np.pi*f*a/C)**4
        self.assertAlmostEqual(abs(sphere(f, 2*a))**2/expected, 1, places=5)
        self.assertAlmostEqual(abs(sphere(79.5e9, .006)-sphere(79.5e9, .006, 12)), 0, places=12)

    def test_instantaneous_frequency(self):
        e = echo(self.p, (.5, 6285, .2), gate=False)
        numerical = np.gradient(e['phase'], 1/self.p['fs'], axis=1)/(2*np.pi)
        self.assertLess(np.max(abs(numerical[:, 2:-2]-e['if_hz'][:, 2:-2])), 10)

    def test_radial_doppler_sign_and_factor_two(self):
        p = dict(self.p, freq_slope=0)
        v = 6700.
        for x0, radial_velocity in [(-1., v), (2., -v)]:
            result = echo(p, (x0, v, 0.), chirps=1, gate=False)
            expected = -2*p['f_start']*radial_velocity/(C+radial_velocity)
            np.testing.assert_allclose(result['if_hz'], expected, rtol=1e-12)
            derivative = np.gradient(result['phase'], 1/p['fs'], axis=1)/(2*np.pi)
            np.testing.assert_allclose(derivative, expected, rtol=1e-8)

    def test_y_sign_unobservable(self):
        a = echo(self.p, (.5, 6285, .2))['iq']
        b = echo(self.p, (.5, 6285, -.2))['iq']
        np.testing.assert_array_equal(a, b)

    def test_sideband_and_finite_path(self):
        approaching = echo(self.p, (2, 6285, .2), chirps=8)
        self.assertEqual(np.max(abs(approaching['iq'])), 0)
        receding = echo(self.p, (-1, 6285, .2))
        self.assertGreater(np.max(abs(receding['iq'])), 0)
        self.assertEqual(np.sum(abs(receding['iq'][16:])**2), 0)

    def test_profiled_amplitude(self):
        h = echo(self.p, (-1, 6285, .2))['iq'].ravel()
        a = 2+3j
        self.assertAlmostEqual(score(a*h, h)/(abs(a)**2*np.vdot(h, h).real), 1, places=12)

    def test_noise_statistic_mean(self):
        rng = np.random.default_rng(25)
        h = np.exp(1j*np.arange(128))
        n = (rng.normal(size=(10000, 128))+1j*rng.normal(size=(10000, 128)))/np.sqrt(2)
        scores = abs(n@h.conj())**2/np.vdot(h, h).real
        self.assertLess(abs(scores.mean()-1), .04)

    def test_complex_search_products(self):
        product = Path(__file__).parent/'search_experiment.h5'
        if not product.exists():
            self.skipTest('Run search_simulated_echoes.py first')
        with h5py.File(product) as h:
            for name in ['parallel', 'perpendicular']:
                g = h[name]
                for ds in ['exact_echo', 'noisy_iq']:
                    z = g[ds][:]
                    self.assertTrue(np.iscomplexobj(z))
                    self.assertGreater(np.max(abs(z.imag)), 0)
                for n in [1,2,4,8]:
                    z = g[f'case_{n}/fitted_echo'][:]
                    self.assertEqual(z.shape, (n,128))
                    self.assertTrue(np.iscomplexobj(z))
                    self.assertGreater(np.max(abs(z.imag)), 0)


if __name__ == '__main__': unittest.main()
