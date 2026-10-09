"""Check browser PEC Mie inversion and radar-equation conversion against SciPy."""
import json
import subprocess
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn, spherical_yn
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
C=299792458.
def reference(frequency,diameter):
    k=2*np.pi*frequency/C
    x=k*diameter/2
    n=np.arange(1,int(np.ceil(x+4*np.cbrt(x)+12))+1)
    j,y=spherical_jn(n,x),spherical_yn(n,x)
    psi=x*j
    xi=x*(j+1j*y)
    dp=j+x*spherical_jn(n,x,derivative=True)
    dx=j+1j*y+x*(spherical_jn(n,x,derivative=True)+1j*spherical_yn(n,x,derivative=True))
    return np.pi/k**2*abs(np.sum((2*n+1)*(-1.)**n*(-dp/dx+psi/xi)))**2

diameters=np.geomspace(.00001,.1,140)
frequency=77e9
sigma=reference(frequency,.002)
parameters=dict(temperature=9000,range=.2,frequency=frequency,time=144e-6,txPowerDbm=12,txGainDbi=6,rxGainDbi=6,lossDb=0)
payload=dict(diameters=diameters.tolist(),frequency=frequency,sigma=float(sigma),parameters=parameters)
script="""import fs from 'node:fs';import {sphereRcs,diameterRoots,estimateRcs} from './web/lab/rcs.mjs';
const p=JSON.parse(fs.readFileSync(0,'utf8'));
console.log(JSON.stringify({rcs:p.diameters.map(d=>sphereRcs(p.frequency,d)),roots:diameterRoots(p.sigma,p.frequency,.00001,.02),single:estimateRcs(101,p.parameters),beam:estimateRcs(401,p.parameters,4)}));"""
result=json.loads(subprocess.check_output(['node','--input-type=module','-e',script],input=json.dumps(payload).encode(),cwd=ROOT))
np.testing.assert_allclose(result['rcs'],[reference(frequency,d) for d in diameters],rtol=2e-10,atol=1e-25)
grid=np.geomspace(.00001,.02,16385)
values=[reference(frequency,d)-sigma for d in grid]
roots=[brentq(lambda d:reference(frequency,d)-sigma,a,b,xtol=1e-15) for a,b,fa,fb in zip(grid[:-1],grid[1:],values[:-1],values[1:]) if fa*fb<0]
np.testing.assert_allclose(result['roots'],roots,rtol=1e-9,atol=1e-12)
assert len(roots)>1,'Test a nonunique Mie inversion'
expected=100*1.380649e-23*9000/144e-6*(4*np.pi)**3*.2**4/(10**((12-30)/10)*10**(12/10)*(C/frequency)**2)
np.testing.assert_allclose([result['single'],result['beam']],expected,rtol=1e-12)
print('PASS: PEC Mie series vs SciPy over 0.01–100 mm, all inversion branches, thermal noise-floor subtraction, and four-RX gain scaling')
