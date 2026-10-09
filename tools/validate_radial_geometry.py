"""Verify the slide geometry identities and monostatic acceleration grid bound."""
import numpy as np
import sympy as sp

t, V, x0, y0 = sp.symbols('t V x0 y0', real=True)
R = sp.sqrt((V*t-x0)**2+y0**2)
assert sp.simplify(sp.diff(R,t)-V*(V*t-x0)/R) == 0
assert sp.simplify(sp.diff(R,t,2)-V**2*y0**2/R**3) == 0

wavelength = 299792458./77e9
epsilon = .2
for H in [10e-6, 100e-6, 1e-3]:
    step = epsilon*wavelength/(2*np.pi*H**2)
    h = np.linspace(-H,H,1001)
    separation = 2*np.pi*step*h**2/wavelength
    assert np.max(separation) <= epsilon*(1+1e-14)
    # Nearest-grid mismatch is at most half an adjacent-trial separation.
    assert np.max(separation/2) <= epsilon/2*(1+1e-14)
    longer_step = epsilon*wavelength/(2*np.pi*(2*H)**2)
    np.testing.assert_allclose(step/longer_step,4.)
print('PASS: exact radial derivatives, endpoint acceleration-phase bound, and quadratic duration scaling')
