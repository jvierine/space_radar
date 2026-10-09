"""Additional regression checks of train-cell guarantees and gap handling."""
import itertools
import numpy as np
from memo_006_train_bank import make_bank, output, phase

f0, slope, fs = 77e9, 99.987e12, 12.5e6
u = 2e-6+np.arange(225)/fs
rng = np.random.default_rng(713)
for bounds,starts in [
    (np.array([[.099,.101],[340,343],[.149,.151]]),np.array([0.,25.37e-6,70e-6])),
    (np.array([[-.101,-.099],[-343,-340],[.149,.151]]),np.array([0.,25.37e-6,70e-6])),
    (np.array([[.099,.101],[-1,1],[-.151,-.149]]),np.array([0.,25.37e-6,70e-6])),
]:
    boxes,centres,errors=make_bank(bounds,starts,u,f0,slope,.112,200000)
    t=starts[:,None]+u; tr=(t.min()+t.max())/2; ur=(u[0]+u[-1])/2
    for i in rng.choice(len(boxes),min(25,len(boxes)),replace=False):
        centre=centres[i];box=boxes[i]
        base=phase(centre,t,u,f0,slope)-phase(centre,tr,ur,f0,slope)
        for corner in itertools.product([0,1],repeat=3):
            p=box[np.arange(3),corner]
            delta=phase(p,t,u,f0,slope)-phase(p,tr,ur,f0,slope)-base
            assert np.max(abs(delta))<=errors[i]+1e-8
    truth=bounds.mean(axis=1); d=np.exp(1j*phase(truth,t,u,f0,slope))
    z=output(d,truth,starts,u,f0,slope,fs,32768)
    assert abs(z/d.size)**2>.999
    print('PASS signed geometry/irregular gaps:',len(boxes),'cells')
# Failure returns no partial bank.
try:
    make_bank(np.array([[.09,.11],[300,400],[.14,.16]]),np.arange(4)*25.37e-6,u,f0,slope,.01,1)
except RuntimeError as e:
    assert 'incomplete' in str(e)
else:
    raise AssertionError('A truncated bank was returned')
try:
    make_bank(np.array([[0,.1],[0,400],[-.1,.1]]),np.arange(4)*25.37e-6,u,f0,slope,.1,100)
except ValueError as e:
    assert 'zero range' in str(e)
else:
    raise AssertionError('Collision geometry was accepted')
print('PASS template-limit and collision failures')
