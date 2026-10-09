"""Compiled-Wasm multi-chirp FFT search, exact scoring, and runtime comparison.

Synthetic inputs only. Run: conda run -n base python tools/validate_fft_search.py
Products: web/lab/qa/fft_search_validation.h5 (not CSV or measurement data).
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess
import time
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("reference", ROOT / "capture_analysis/analyze_captures.py")
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)
p = dict(rows=1000, n=225, frame=125, fs=12.5e6, adc=2e-6,
         period=25.37e-6, f0=76999999970.19768, slope=99987387657165.53)
profile = dict(samples=p["n"], T_adc=p["adc"], period=p["period"], fs=p["fs"],
               f_start=p["f0"], freq_slope=p["slope"], if_max=10e6)


def template(theta, pulses):
    x, y, v = theta
    echo = reference.echo(profile, [x, v, y], pulses, gate=False, phase_only=True)
    return (np.exp(1j * echo["phase"]) / echo["range"] ** 2 *
            (echo["u"] >= 2 * echo["range"] / reference.C))


def correction(theta, pulses):
    x,y,v=theta
    u=p["adc"]+np.arange(p["n"])/p["fs"]
    m=(u[0]+u[-1])/2
    d=(np.arange(pulses)-(pulses-1)/2)*p["period"]
    star=m+(pulses-1)*p["period"]/2
    beta=v/reference.C
    X=x-v*star
    r=(np.sqrt(X*X+(1-beta*beta)*y*y)+beta*X)/(1-beta*beta)
    rd=-v*(X/np.sqrt(X*X+(1-beta*beta)*y*y)+beta)/(1-beta*beta)
    slow=-2*(p["f0"]+p["slope"]*m)*rd/reference.C+4*p["slope"]*r*rd/reference.C**2
    fast=slow-2*p["slope"]*r/reference.C
    return template(theta,pulses)*np.exp(-2j*np.pi*(fast*(u-m)+slow*d[:,None]))


rng = np.random.default_rng(63)
sigma = 20.0
z = (rng.normal(size=(1000, 225)) + 1j * rng.normal(size=(1000, 225))) * sigma / np.sqrt(2)
z += 400 * np.exp(2j * np.pi * np.arange(225) / 23)
truth = template([.3, .1, 341.51], 16)
z[750:766] += truth / np.linalg.norm(truth) * sigma * np.sqrt(3000)
z[880:884] = 0
z = z.astype(np.complex64)
output = ROOT / "web/lab/qa/fft_search_validation.h5"
output.parent.mkdir(parents=True, exist_ok=True)
wasm_hash = hashlib.sha256((ROOT / "web/lab/core.wasm").read_bytes()).hexdigest()
with h5py.File(output, "w") as h:
    h.attrs["wasm_sha256"] = wasm_hash
    h.attrs["generator"] = "tools/validate_fft_search.py + tools/benchmark_fft.mjs"
    h.attrs["input"] = "Synthetic retarded Eq.25 echo plus complex Gaussian noise and static background"
    h.attrs["limitations"] = "Timing benchmarks are machine-specific. Local sampled grid accuracy is not a global certificate. Synthetic recovery is not a measured detection."
    for name, data, bounds, pulses in [
        ("focused", z, [.29, .31, .09, .11, 331.51, 351.51], [4, 8, 16]),
        ("broad", z, [-.6, .6, .05, .3, 200, 500], [4]),
        ("mirrored_iq", z.conj(), [-.31, -.29, .09, .11, -351.51, -331.51], [4, 8, 16]),
    ]:
        started = time.perf_counter()
        run = subprocess.run(["node", "tools/benchmark_fft.mjs"],
            input=json.dumps(dict(p=p, data=np.stack([data.real, data.imag], -1).ravel().tolist(),
                                  bounds=bounds, pulses=pulses)), text=True, capture_output=True, check=True, cwd=ROOT)
        print(run.stderr, end="", flush=True)
        group = h.create_group(name)
        group.attrs["bounds_x_y_v"] = np.asarray(bounds).reshape(3, 2)
        group.attrs["harness_wall_seconds"] = time.perf_counter() - started
        for row in json.loads(run.stdout):
            best, direct, verified = map(np.asarray, [row["best"], row["direct"], row["verified"]])
            np.testing.assert_allclose(best[3], verified[3], rtol=5e-5)
            assert np.all(best[:3] >= np.array(bounds)[::2] - 1e-5)
            assert np.all(best[:3] <= np.array(bounds)[1::2] + 1e-5)
            assert int(best[4]) == (1 if name == "mirrored_iq" else 0)
            n = row["pulses"]
            q = template(best[:3], n).ravel()
            expected_truth = truth[:n].ravel()
            coherence = abs(np.vdot(q, expected_truth))**2 / (np.vdot(q,q).real * np.vdot(expected_truth,expected_truth).real)
            direct_q = template(direct[:3], n).ravel()
            direct_coherence = abs(np.vdot(direct_q,expected_truth))**2/(np.vdot(direct_q,direct_q).real*np.vdot(expected_truth,expected_truth).real)
            assert coherence > .9, (name,n,coherence,direct_coherence)
            initial_nodes = np.asarray(row["initial_nodes"]).reshape(-1,3)
            if name != "broad":
                truth_correction=correction([.3,.1,341.51],n).ravel()
                coverage=0.
                for th in initial_nodes:
                    c=correction(th,n).ravel()
                    coverage=max(coverage,abs(np.vdot(c,truth_correction))**2/(np.vdot(c,c).real*np.vdot(truth_correction,truth_correction).real))
                assert coverage > .94, (name,n,coverage)
            # Do not require exact coordinates: short trains have trajectory ridges.
            assert 10*np.log10(best[3]/direct[3]) > -.2, (name,n,best,direct)
            values = np.asarray(row["fitted"]).reshape(-1, 2)
            fit = values[:, 0] + 1j * values[:, 1]
            if best[4]: q = q.conj()
            residual = data.astype(np.complex128) - data[:125].astype(np.complex128).mean(0)
            observed = residual[750:750+n].ravel()
            expected_fit = q * np.vdot(q, observed) / np.vdot(q, q)
            np.testing.assert_allclose(fit, expected_fit, rtol=2e-4, atol=.002)
            nodes = np.asarray(row["nodes"]).reshape(-1, 3)
            np.testing.assert_allclose(np.nanmax(row["scores"]), best[3], rtol=1e-6)
            assert np.all(nodes >= np.asarray(bounds)[::2] - 1e-5)
            assert np.all(nodes <= np.asarray(bounds)[1::2] + 1e-5)
            saved = group.create_group(str(n))
            saved.attrs["template_coherence"] = coherence
            if name != "broad": saved.attrs["correction_coverage_coherence"] = coverage
            saved.attrs["fft_total_seconds"] = sum(row[k] for k in ["build_seconds", "fft_seconds", "refine_seconds"])
            saved.attrs["direct_seconds"] = row["direct_seconds"]
            for key in ["specification", "best", "direct", "verified", "boxes", "nodes", "initial_nodes", "scores", "fitted"]:
                saved.create_dataset(key, data=np.asarray(row[key]), compression="gzip")
            print(f"{name}: {n} chirps; coherence {coherence:.6f}; FFT/direct peak difference {10*np.log10(best[3]/direct[3]):+.3f} dB", flush=True)
print(f"PASS: multi-chirp searches and exact scores/fits; {output}")
