"""Reproducible FYS-3000 space-radar numbers and Eq.25 teaching records.

Run with conda run -n base python space_radar_assets.py. No CSV products.
Uses the existing conducting-sphere solver and the exact memo simulation.
"""
from pathlib import Path
import importlib.util
import hashlib
import json

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
CAPTURE = HERE.parent / 'capture_analysis'
spec = importlib.util.spec_from_file_location('lecture_echo', CAPTURE / 'analyze_captures.py')
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)
C, KB = model.C, model.KB
PT, TSYS = 10**((12-30)/10), 10000.
N = np.array([1, 2, 4, 10])
R = np.array([1., 10., 100.])


def main():
    out = HERE / 'assets/fmcw_lecture.h5'
    with h5py.File(CAPTURE / 'search_experiment.h5', 'r') as old, h5py.File(out, 'w') as h:
        h.attrs['generator'] = Path(__file__).name
        h.attrs['generator_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        h.attrs['model_path'] = str(CAPTURE / 'analyze_captures.py')
        h.attrs['model_sha256'] = hashlib.sha256((CAPTURE / 'analyze_captures.py').read_bytes()).hexdigest()
        h.attrs['memo_signal_sha256'] = hashlib.sha256((CAPTURE / 'search_experiment.h5').read_bytes()).hexdigest()
        h.attrs['assumptions'] = 'Maximum power confirmed by user: 12 dBm per TX; one TX and one RX; perfectly conducting sphere; Tsys=10000 K; gain sensitivity 6 and 10 dBi per antenna; extra loss 0 dB'
        x = np.geomspace(.001, 10., 700)
        # Set wavelength=1 m: x denotes DIAMETER / wavelength, not radius.
        sigma = abs(model.sphere(C, x))**2
        area = np.pi*(x/2)**2
        normalized = sigma/area
        rayleigh = 9*(np.pi*x)**4
        assert abs(normalized[0]/rayleigh[0]-1) < 1e-3
        g = h.create_group('scattering')
        g.create_dataset('diameter_over_wavelength', data=x)
        g.create_dataset('mie_normalized_rcs', data=normalized)
        g.create_dataset('rayleigh_normalized_rcs', data=rayleigh)
        freqs = np.array([3e9, 10e9, 78e9, 79.5e9])
        sigma3 = abs(model.sphere(freqs, .003))**2
        g.create_dataset('frequency_hz', data=freqs)
        g.create_dataset('three_mm_rcs_m2', data=sigma3)

        lecture = h.create_group('three_mm')
        lecture.create_dataset('chirps', data=N)
        lecture.create_dataset('range_m', data=R)
        fs, samples, period = 12.5e6, 128, 13.709999904632568e-6
        lecture.attrs['frequency_hz'] = 78e9
        lecture.attrs['rcs_m2'] = float(sigma3[2])
        lecture.attrs['tx_power_w'] = PT
        lecture.attrs['system_temperature_K'] = TSYS
        lecture.attrs['adc_on_time_per_chirp_s'] = samples/fs
        lecture.attrs['assumed_chirp_period_s'] = period
        live = N*samples/fs
        lecture.create_dataset('adc_on_time_s', data=live)
        lecture.create_dataset('window_span_s', data=(N-1)*period+samples/fs)
        lecture.create_dataset('equivalent_noise_bandwidth_hz', data=1/live)
        for gain in [6, 10]:
            G = 10**(gain/10)
            power = PT*G*G*(C/78e9)**2*sigma3[2]/((4*np.pi)**3*R**4)
            snr = power[None,:]*live[:,None]/(KB*TSYS)
            lecture.create_dataset(f'snr_gain_{gain}_db', data=10*np.log10(snr))
        antenna = h.create_group('antenna')
        gain = np.array([6., 10.])
        G = 10**(gain/10)
        omega = 4*np.pi/G
        cone = 2*np.rad2deg(np.arccos(1-omega/(2*np.pi)))
        antenna.create_dataset('gain_dbi', data=gain)
        antenna.create_dataset('linear_gain', data=G)
        antenna.create_dataset('ideal_lossless_beam_solid_angle_sr', data=omega)
        antenna.create_dataset('ideal_uniform_cone_full_angle_deg', data=cone)
        antenna.create_dataset('effective_aperture_mm2', data=G*(C/78e9)**2/(4*np.pi)*1e6)

        for name in ['parallel', 'perpendicular']:
            p = json.loads(old[name].attrs['profile_json'])
            truth = old[name]['truth'][:]
            e = model.echo(p, truth, 32)
            np.testing.assert_allclose(e['iq'], old[name]['exact_echo'][:], rtol=1e-12, atol=1e-20)
            group = h.create_group(name)
            group.attrs['profile_json'] = json.dumps(p)
            group.attrs['geometry_status'] = 'Same illustrative coordinates as memo; not a chamber survey'
            group.attrs['maximum_power_confirmed'] = True
            group.create_dataset('truth', data=truth)
            for key,val in e.items(): group.create_dataset(key, data=val)
            sig = abs(model.sphere(p['f_start'], p['projectile_size']*.001))**2
            power = PT*10**(12/10)*(C/p['f_start'])**2*sig/(4*np.pi)**3
            group.attrs['rcs_m2'] = sig
            group.attrs['one_rx_power_1m_W_gain6'] = power
            group.create_dataset('ideal_1m_snr_db',data=10*np.log10(power*N*samples/fs/(KB*TSYS)))
            energy = np.sum(abs(e['iq'])**2,axis=1)/fs
            moving = np.array([np.max(np.convolve(energy,np.ones(n),'valid'))/(KB*TSYS) for n in N])
            group.create_dataset('illustrative_moving_best_window_snr_db',data=10*np.log10(moving))
            group.attrs['assumed_1m_visibility_s'] = 1/p['speed']
            group.attrs['assumed_1m_chirp_periods'] = 1/(p['speed']*p['period'])
            complete_intervals=(1/p['speed']-samples/fs)/p['period']
            group.attrs['assumed_1m_full_adc_chirps_min'] = int(np.floor(complete_intervals))
            group.attrs['assumed_1m_full_adc_chirps_max'] = int(np.ceil(complete_intervals))
            group.attrs['visibility_count_convention'] = 'Dwell/chirp period is an approximate count; fully contained ADC records depend on arbitrary entry phase. Assumes continuous chirp train over 1m of useful trajectory, no frame gap.'

        p = json.loads(old['parallel'].attrs['profile_json'])
        truth = old['parallel/truth'][:]
        signal = old['parallel/exact_echo'][:3]
        velocities = np.array([6200., 7200., truth[1]])
        candidates, products, integrals, coherence = [], [], [], []
        for velocity in velocities:
            theta=truth.copy();theta[1]=velocity
            q=model.echo(p,theta,3)['iq']
            candidates.append(q)
            product = signal*q.conj()
            products.append(product)
            integrals.append(np.cumsum(product.ravel())/p['fs'])
            coherence.append(abs(np.vdot(q,signal))**2/(np.vdot(q,q).real*np.vdot(signal,signal).real))
        assert np.isclose(coherence[-1], 1., atol=1e-12)
        assert max(coherence[:-1]) < .5
        mf = h.create_group('matched_filter_demo')
        mf.attrs['source'] = 'First three consecutive chirps of parallel/exact_echo in memo search_experiment.h5'
        mf.create_dataset('signal_iq',data=signal)
        mf.create_dataset('time_s',data=model.echo(p,truth,3)['t'])
        mf.create_dataset('candidate_velocity_m_s',data=velocities)
        mf.create_dataset('template_iq',data=candidates)
        mf.create_dataset('mixed_iq',data=products)
        mf.create_dataset('cumulative_integral',data=integrals)
        mf.create_dataset('power_coherence',data=coherence)
        mf.attrs['signal_energy_J'] = float(np.sum(abs(signal)**2)/p['fs'])
        mf.attrs['adc_on_time_s'] = 3*samples/fs
        mf.attrs['window_span_s'] = 2*p['period']+samples/fs
        mf.attrs['equivalent_noise_bandwidth_hz_constant_amplitude'] = fs/(3*samples)
        print('3 mm RCS:',dict(zip(freqs/1e9,sigma3)),flush=True)
        print('3 mm, gains 6 dBi, rows 1/2/4/10 chirps, columns 1/10/100 m:',lecture['snr_gain_6_db'][:],flush=True)
        print('MF candidate coherence:',coherence,flush=True)
        for name in ['parallel','perpendicular']:
            print(name,'ideal 1m SNR',h[name]['ideal_1m_snr_db'][:], 'moving illustration',h[name]['illustrative_moving_best_window_snr_db'][:],flush=True)
    print(out,flush=True)


if __name__ == '__main__':
    main()
