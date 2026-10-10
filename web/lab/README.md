# Coherent FMCW data analysis

Live application: **[juha.no/fmcw/lab/](https://juha.no/fmcw/lab/)**.
The [repository README](../../README.md#1-coherent-integration-data-analysis-rust--wasm--webgpu)
maps every main code component and gives build/test commands.

## Use

1. Move/resize the blue **Background** window to a quiet interval.
2. Move/resize the purple **Analysis** window.
3. Choose coherent chirp count and search bounds; press **Play**.
4. Watch results accumulate. **Download HDF5** enables after a completed scan.

Drag empty colormesh space to zoom; **Reset full view** restores the recording.
Moving the yellow coherent train runs that selected train's fit. Shared `?gui=...`
URLs preserve settings and view. The bare URL loads the configured default
windows, eight chirps and r0/v0/a0 bounds. Voltage and spectrum use the selected
receiver; trajectory matched powers use all four receivers.

## Build and data import

Run from the repository root:

```sh
rustup target add wasm32-unknown-unknown
./tools/build_lab.sh
conda run -n base python tools/prepare_lab_data.py /path/to/recording.nc
conda run --no-capture-output -n base python -m http.server 18791 --bind 127.0.0.1 --directory web
```

Open http://127.0.0.1:18791/lab/. The compiled Wasm is tracked; rebuilding is
necessary only for Rust changes. Recording files are not bundled. The importer
needs NumPy, h5py and netCDF4 and writes HDF5 plus lossless interleaved float32
transport under `web/lab/datasets/test63/`, ignored by Git. Original deployed
data live under `/mnt/shovel/fmcw/` on `juha-no`.

## Processing

- Rust/Wasm owns the automatic r0/v0/a0 grid, templates, sparse FFT bank and
  exact candidate refinement. WebGPU runs complex64 FFT/scoring kernels, with
  numerical checks and CPU fallback. Analysis runs in a worker, in the browser.
- `R(t)=r0+v0*h+0.5*a0*h²`, `h=t-t_mid`. Positive v0 means increasing slant
  range. a0 is radial range curvature; the polynomial is a local approximation.
  Acceleration planes and bounded velocity groups share phase corrections;
  fast-time and chirp-time FFTs search the dominant sinusoids.
- Each RX has its own complex quiet mean at each fast-time sample. Noise comes
  from the **mean-removed quiet voltage**, using all intact background samples
  and correcting variance for fitting that mean. Stationary wall power is removed.
  One background chirp only permits a raw-power upper bound; RCS is then omitted.
- Full-bandwidth RX noise is scaled to `B_analysis=1/T_coh`, where
  `T_coh=chirps*samples/fs`. Unsampled gaps enter phase timing, not acquired energy.
  Excess-signal SNR displays with a 0 dB floor; raw ratios and physical powers
  remain unchanged. See [GPU.md](GPU.md) for normalization details.
- The trajectory bank sums RX matched powers incoherently and searches both
  I/Q orientations. Maps use **MAX** over the omitted parameter. The 32-million
  grid-node cap rejects excessive banks rather than narrowing supplied bounds.
- Antenna phasing uses a relative-phase grid and Nelder–Mead refinement. Three
  phase projections, individual/beam SNR, RCS and all sphere-diameter roots are
  displayed. See [BEAMFORMING.md](BEAMFORMING.md) for assumptions.
- Waveforms and fitted continuous range are lines. Analyzed histories are
  unconnected scatter points. The HDF5 export includes completed train results,
  phase/SNR/RCS/diameter estimates and the current full bank and voltage fit.
- Padding/nonfinite chirps are excluded; coherent trains remain within a frame.
  Times are reconstructed seconds since file start. Neither a fit nor a search
  maximum alone confirms a projectile detection.

## Validation

See [GPU tests, fallback and benchmark provenance](GPU.md). Additional regressions:

```sh
node tools/validate_background_residual_noise.mjs
node tools/validate_selected_spectrum.mjs
node tools/validate_analysis_controls.mjs
node tools/validate_scan_history.mjs
node tools/validate_analysis_export.mjs
node tools/validate_snr_display.mjs
node tools/validate_complexity_counts.mjs
```

Commands run from the repository root; measurement-dependent tests need prepared
Test 63 data. Python commands use base conda; scientific outputs use HDF5.

## Coherent time–range map

Playback adds one range column per completed coherent train. The best joint
trajectory fit supplies velocity, acceleration, and four complex receive beam
weights; all three remain fixed across the column. Each range is evaluated
with the full FMCW template, using complex antenna summation and temporal
coherent integration. Noise is propagated as `(qᴴq)(wᴴCw)` from the quiet
receive covariance. Color is excess matched SNR in dB, floored at 0 dB;
selection of the peak and its weights can bias noisy estimates. Gaps remain
blank. The current range profile and all profile powers, noise powers, scores,
weights and motion parameters are included in the HDF5 download.

## Recording selector

The recording dropdown reads `dataset-catalog.json`, with labels such as
`Shot 63 · 2026-07-23`. Each entry has an ID and a `datasets/<id>/` path.
Switching records reloads the workers and clears the old analysis and plot
state. Shared links retain the selected recording through `?dataset=<id>`.
Sampling rate, chirp dimensions and timing come from that recording's metadata.

Package another NetCDF recording using:

```sh
conda run -n base python tools/prepare_lab_data.py /path/to/shot.nc --output web/lab/datasets/shot99 --id shot99 --shot 99
```

Add its ID, shot/date label and path to the catalog only after packaging and
verifying the recording. Explicit frame timing can be supplied with
`--frame-interval-s`; `--ramp-tail-us` records the reconstruction assumption.
Recording arrays remain outside Git.

The hypervelocity records are Shot 6606 (2025-02-24) and Shot 6607
(2025-02-25), supplied by Pål Ellingsen. Browser inputs retain original
complex ADC counts for source frames 43–46 and 34–37 respectively, with
full-source SHA-256 and source frame offsets recorded in metadata/HDF5.
The header identifies these event excerpts. Reconstructed times retain the
source frame offset; frame gaps remain unknown unless explicit timing is
supplied. Velocity defaults span −7 to +7 km/s; no search runs automatically.
The original recordings remain on the Expansion drive.

For an excerpt, add `--frame-start 43 --frame-stop 47` to the importer command;
frame-stop is exclusive. The importer supports both Test 63 attributes and
the earlier hypervelocity files' `projectile_size` variable.

After importing both excerpts, run `conda run -n base python
tools/configure_hypervelocity_datasets.py` to apply the shot-specific quiet
windows, display selections and velocity bounds and register both records.
