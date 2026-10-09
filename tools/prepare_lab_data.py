"""Package original NetCDF as HDF5 plus lossless browser transport arrays.
Run in base conda. The .f32 files are HTTP/Wasm transport, not archival products.
"""

import argparse, hashlib, json
from pathlib import Path
import h5py, numpy as np
from netCDF4 import Dataset, num2date

p = argparse.ArgumentParser()
p.add_argument("source", type=Path)
p.add_argument("--output", type=Path, default=Path("web/lab/datasets/test63"))
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
sha = hashlib.sha256(a.source.read_bytes()).hexdigest()
with Dataset(a.source, auto_complex=True) as source:
    raw = np.asarray(source["radar_cube"][0], dtype=np.complex64)
    frames, chirps, rx, samples = raw.shape
    scalar = {
        key: float(source[key][:])
        for key in ["fs", "freq_slope", "f_start", "T_adc", "T_idle", "if_max", "speed"]
    }
    t = source["time"]
    clock = num2date(t[0], t.units, calendar=t.calendar).strftime("%Y-%m-%dT%H:%M:%S")
    meta = dict(
        id="test63",
        filename=a.source.name,
        source_sha256=sha,
        frames=frames,
        chirps_per_frame=chirps,
        receivers=rx,
        samples=samples,
        total_chirps=frames * chirps,
        recorded_clock=clock,
        clock_timezone="not specified in source",
        frame_of_interest=int(source["frame_of_interest"][:]),
        parameters=scalar,
        description=str(source.description),
        orientation=str(source.orientation),
        material=str(source.material),
        diameter_mm=float(source.ball_size_mm),
        zero_padded_samples=int(source.num_zero_padded_samples),
        ramp_tail_us_assumption=0.37,
        frame_interval_s=None,
        transport=[],
    )
    with h5py.File(a.output / "recording.h5", "w") as archive:
        archive.create_dataset(
            "complex_counts", data=raw, compression="gzip", shuffle=True
        )
        archive.attrs["source_sha256"] = sha
        archive.attrs["source_filename"] = a.source.name
        archive.attrs["generator"] = "tools/prepare_lab_data.py"
        archive.attrs["generator_sha256"] = hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest()
        archive.attrs["recorded_clock"] = clock
        archive.attrs["clock_timezone"] = meta["clock_timezone"]
        archive.attrs["parameters_json"] = json.dumps(scalar)
    for k in range(rx):
        z = raw[:, :, k, :].reshape(frames * chirps, samples).copy()
        stream = np.stack([z.real, z.imag], axis=-1).astype("<f4")
        path = a.output / f"rx{k}.f32"
        path.write_bytes(stream.tobytes())
        assert np.array_equal(
            np.frombuffer(path.read_bytes(), dtype="<f4").reshape(*z.shape, 2)[..., 0]
            + 1j
            * np.frombuffer(path.read_bytes(), dtype="<f4").reshape(*z.shape, 2)[
                ..., 1
            ],
            z,
        )
        meta["transport"].append(
            dict(
                rx=k,
                file=path.name,
                bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            )
        )
    (a.output / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
print(json.dumps(meta, indent=2))
