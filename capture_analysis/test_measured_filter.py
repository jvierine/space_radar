"""Reference checks of the coloured-noise, per-RX measured-data statistic."""
import unittest
import numpy as np
from apply_measured_filter import MeasuredBank, prepare, rail_mask, windows, NC


class MeasuredStatisticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=dict(fs=12.5e6,freq_slope=60e12,f_start=79.5e9,T_adc=0,
                   if_max=10e6,period=13.71e-6,samples=8,projectile_size=6,speed=6285)
        cls.good=np.ones((2,8),bool)
        cls.good[:,0]=False
        cls.good[1,3]=False
        cls.psd=np.array([np.linspace(1,3,8),np.linspace(4,2,8)])
        cls.bank=MeasuredBank(cls.p,cls.good,cls.psd)

    def test_coarse_matches_direct_projection_with_mask_and_colour(self):
        theta=np.array([-.1,6285,.2])
        h=self.bank.base.template(theta).reshape(NC,8)
        gains=np.array([2+3j,-1+.5j])
        counts=h[:,None,:]*gains[None,:,None]
        white=prepare(counts,np.zeros((2,8)),self.good,self.psd)
        direct,fitted,recovered=self.bank.evaluate(theta,white,0,full=True)
        np.testing.assert_allclose(recovered,gains,rtol=1e-12)
        self.assertAlmostEqual(direct/np.sum(abs(white)**2),1,places=12)
        np.testing.assert_allclose(fitted,counts,rtol=1e-12)
        coarse=self.bank.coarse(white,np.array([0]))
        idx=int(np.argmin(np.linalg.norm((self.bank.parameters-theta)/self.bank.step,axis=1)))
        self.assertLess(abs(coarse[0,idx,0]/direct-1),2e-6)
        # A reference must check which sampled I/Q orientation is retained.
        flipped=prepare(counts.conj(),np.zeros((2,8)),self.good,self.psd)
        self.assertAlmostEqual(self.bank.evaluate(theta,flipped,1)/np.sum(abs(flipped)**2),1,places=12)

    def test_rail_flag_and_selection_stays_within_frame(self):
        z=np.array([1+2j,-32768+3j,4+32767j])
        np.testing.assert_array_equal(rail_mask(z),[False,True,True])
        rng=np.random.default_rng(4)
        white=rng.normal(size=(1030,2,8))+1j*rng.normal(size=(1030,2,8))
        starts,energy,onset=windows(white,self.p['fs'])
        self.assertTrue(np.all(starts>=0))
        self.assertTrue(np.all(starts+NC<=1030))
        self.assertEqual(len(energy),1030)


if __name__=='__main__': unittest.main()
