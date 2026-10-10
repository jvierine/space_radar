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
p.add_argument("--id", help="Catalog recording ID; defaults to output folder name")
p.add_argument("--shot", help="Shot number for the recording label")
p.add_argument("--ramp-tail-us", type=float, default=0.37)
p.add_argument("--frame-interval-s", type=float)
p.add_argument("--frame-start", type=int, default=0)
p.add_argument("--frame-stop", type=int)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8*1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()
sha = digest(a.source)
with Dataset(a.source, auto_complex=True) as source:
    source_frames = source["radar_cube"].shape[1]
    stop = a.frame_stop if a.frame_stop is not None else source_frames
    if not 0 <= a.frame_start < stop <= source_frames:
        raise ValueError("Invalid source frame interval")
    raw = np.asarray(source["radar_cube"][0, a.frame_start:stop], dtype=np.complex64)
    frames, chirps, rx, samples = raw.shape
    scalar = {
        key: float(source[key][:])
        for key in ["fs", "freq_slope", "f_start", "T_adc", "T_idle", "if_max", "speed"]
    }
    t = source["time"]
    clock = num2date(t[0], t.units, calendar=t.calendar).strftime("%Y-%m-%dT%H:%M:%S")
    meta = dict(
        id=a.id or a.output.name,
        shot_number=a.shot,
        label=f"Shot {a.shot} · {clock[:10]}" if a.shot else clock[:10],
        filename=a.source.name,
        source_sha256=sha,
        frames=frames,
        chirps_per_frame=chirps,
        receivers=rx,
        samples=samples,
        total_chirps=frames * chirps,
        recorded_clock=clock,
        clock_timezone="not specified in source",
        frame_of_interest=max(0, int(source["frame_of_interest"][:])-a.frame_start),
        source_frame_start=a.frame_start,
        source_frame_stop=stop,
        source_total_frames=source_frames,
        parameters=scalar,
        description=str(getattr(source,"description",getattr(source,"summary",""))),
        orientation=str(getattr(source,"orientation",getattr(source,"comment",""))),
        material=str(getattr(source,"material","not specified")),
        diameter_mm=float(source.ball_size_mm) if "ball_size_mm" in source.ncattrs() else float(source["projectile_size"][:]),
        zero_padded_samples=int(getattr(source,"num_zero_padded_samples",0)),
        ramp_tail_us_assumption=a.ramp_tail_us,
        frame_interval_s=a.frame_interval_s,
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
        archive.attrs["source_frame_start"] = a.frame_start
        archive.attrs["source_frame_stop"] = stop
        archive.attrs["source_total_frames"] = source_frames
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
