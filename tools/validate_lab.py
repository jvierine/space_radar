"""Independent NumPy/memo comparison for compiled Rust/Wasm DSP. Products: HDF5."""

from pathlib import Path
import importlib.util, json, subprocess, hashlib
import numpy as np, h5py

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "memo", ROOT / "capture_analysis/analyze_captures.py"
)
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)
meta = json.loads((ROOT / "web/lab/datasets/test63/metadata.json").read_text())
par = meta["parameters"]
p = dict(
    rows=1000,
    n=225,
    frame=125,
    fs=par["fs"],
    adc=par["T_adc"],
    period=25.37e-6,
    f0=par["f_start"],
    slope=par["freq_slope"],
)
profile = dict(
    samples=p["n"],
    T_adc=p["adc"],
    period=p["period"],
    fs=p["fs"],
    f_start=p["f0"],
    freq_slope=p["slope"],
    if_max=10e6,
)


def q(x, y, v, n, receiver=False):
    e = model.echo(profile, [x, v, y], n, gate=False, phase_only=True)
    z = np.exp(1j * e["phase"]) / e["range"] ** 2 * (e["u"] >= 2 * e["range"] / model.C)
    if receiver:
        f = e["if_hz"]
        z *= 1j * f / (175e3 + 1j * f) * 1j * f / (350e3 + 1j * f)
    return z


rng = np.random.default_rng(63)
sigma = 20.0
noise = (
    (rng.normal(size=(1000, 225)) + 1j * rng.normal(size=(1000, 225)))
    * sigma
    / np.sqrt(2)
)
bg = 400 * np.exp(2j * np.pi * np.arange(225) / 23)
z = noise + bg
truth = q(0.3, 0.1, 341.51, 16)
z[750:766] += truth / np.linalg.norm(truth) * sigma * np.sqrt(3000)
z[880:884] = 0
z = z.astype(np.complex64)
cases = [
    dict(x=x, y=0.1, v=v, pulses=n, receiver=receiver)
    for x, v in [(0.3, 341.51), (-0.3, -341.51), (0.3, -7000), (0.1, 0)]
    for n in [1, 2, 4, 8, 16]
    for receiver in [False, True]
]
request = dict(p=p, data=np.stack([z.real, z.imag], -1).ravel().tolist(), cases=cases)
run = subprocess.run(
    ["node", str(ROOT / "tools/validate_lab.mjs")],
    input=json.dumps(request),
    text=True,
    cwd=ROOT,
    check=True,
    capture_output=True,
)
actual = json.loads(run.stdout)
mean = z[:125].astype(np.complex128).mean(0)
res = z.astype(np.complex128) - mean
assert actual["baselineCount"] == 125
for case, value in zip(cases, actual["templates"]):
    got = np.asarray(value).reshape(-1, 2)
    got = got[:, 0] + 1j * got[:, 1]
    expected = q(
        case["x"], case["y"], case["v"], case["pulses"], case["receiver"]
    ).ravel()
    np.testing.assert_allclose(got, expected, rtol=2e-6, atol=2e-5)
trace = np.asarray(actual["trace"])
np.testing.assert_array_equal(trace[:225], z[750].real)
np.testing.assert_array_equal(trace[225:450], z[750].imag)
np.testing.assert_allclose(trace[450:675], res[750].real, rtol=2e-6, atol=2e-5)
np.testing.assert_allclose(trace[675:900], res[750].imag, rtol=2e-6, atol=2e-5)
w = np.hanning(225)
rawpower = abs(np.fft.fftshift(np.fft.fft(z[750] * w, 256))) ** 2 / w.sum() ** 2
subpower = abs(np.fft.fftshift(np.fft.fft(res[750] * w, 256))) ** 2 / w.sum() ** 2
np.testing.assert_allclose(trace[900:1156], rawpower, rtol=2e-6, atol=2e-5)
np.testing.assert_allclose(trace[1156:1412], subpower, rtol=2e-6, atol=2e-5)
score_errors = []
for row in actual["scores"]:
    n = row["n"]
    starts = []
    k = 125
    while k + n <= 625:
        if k // 125 == (k + n - 1) // 125:
            starts.append(k)
            k += n
        else:
            k += 1
    count = min(len(starts), 24)
    starts = [starts[i * len(starts) // count] for i in range(count)]
    controls = np.stack([res[s : s + n].ravel() for s in starts])
    expected = []
    for x in np.linspace(0.29, 0.31, 3):
        for y in np.linspace(0.09, 0.11, 3):
            for v in np.linspace(331.51, 351.51, 3):
                a = q(x, y, v, n).ravel()
                scores = []
                for a in [a, a.conj()]:
                    scores.append(
                        abs(np.vdot(a, res[750 : 750 + n].ravel())) ** 2
                        / np.mean(abs(controls @ a.conj()) ** 2)
                    )
                expected.append(max(scores))
    got = np.asarray(row["values"][:27])
    np.testing.assert_allclose(got, expected, rtol=3e-5, atol=1e-5)
    score_errors.append(float(np.max(abs(got - expected))))
    best = int(np.argmax(got))
    ix, iy, iv = np.unravel_index(best, (3, 3, 3))
    fitted = q(
        np.linspace(0.29, 0.31, 3)[ix],
        np.linspace(0.09, 0.11, 3)[iy],
        np.linspace(331.51, 351.51, 3)[iv],
        n,
    ).ravel()
    injected = truth[:n].ravel()
    coherence = abs(np.vdot(fitted, injected)) ** 2 / (
        np.vdot(fitted, fitted).real * np.vdot(injected, injected).real
    )
    print(
        "Injected echo",
        n,
        "pulses: peak",
        10 * np.log10(got.max()),
        "dB; template coherence",
        coherence,
    )
    assert coherence > 0.90, ("poor injected-waveform match", n, coherence)
    ideal_energy_ratio = (
        3000
        * np.vdot(truth[:n].ravel(), truth[:n].ravel()).real
        / np.vdot(truth.ravel(), truth.ravel()).real
    )
    assert got[13] > 0.35 * ideal_energy_ratio, (
        "injected signal energy not recovered",
        n,
        got[13],
        ideal_energy_ratio,
    )
assert actual["paddingRejected"] == -1 and actual["crossFrameRejected"] == -1
scan = np.asarray(actual["scanned"]).reshape(-1, 5)
assert not np.any(scan[:, 0] < 625)
assert not np.any((scan[:, 0] <= 883) & (scan[:, 0] + 4 > 880))
assert abs(int(scan[np.argmax(scan[:, 1]), 0]) - 750) < 4
out = ROOT / "web/lab/qa"
out.mkdir(exist_ok=True)
with h5py.File(out / "wasm_validation.h5", "w") as h:
    h.attrs["status"] = (
        "PASS: 40 retarded templates, complex mean, Hann FFT, 5 pulse-count filter cubes, injected waveform/time recovery, padding and frame boundaries"
    )
    h.attrs["wasm_sha256"] = hashlib.sha256(
        (ROOT / "web/lab/core.wasm").read_bytes()
    ).hexdigest()
    h.attrs["generator"] = "tools/validate_lab.py"
    h.create_dataset("maximum_score_absolute_errors", data=score_errors)
    h.create_dataset("time_scan", data=scan)
    for row in actual["scores"]:
        h.create_dataset(f"score_{row['n']}_pulses", data=row["values"])
print(
    "PASS: 40 templates; complex mean; raw/subtracted complex Hann FFT; 1/2/4/8/16 pulse cubes match NumPy; injected waveform and time recovered; padding/frame boundaries excluded."
)
print(out / "wasm_validation.h5")
