# FMCW radar: analysis, simulator and slides

By **Juha Vierinen**.

This repository contains three applications. **Start with the laboratory for
coherent integration of recorded four-receiver data.**

| Part | Live application | Code |
| --- | --- | --- |
| **Coherent-integration data analysis** | [FMCW laboratory](https://juha.no/fmcw/lab/) | [`lab-core/src/`](lab-core/src/) — Rust DSP; [`web/lab/`](web/lab/) — browser UI and WebGPU |
| Measurement simulator | [Signal simulator](https://juha.no/fmcw/sim/) | [`web/sim/`](web/sim/) |
| Slides | [Radar lecture](https://juha.no/fmcw/slides/), [radial FFT search and GUI](https://juha.no/fmcw/slides/radial-fft/) | [`manim/`](manim/) |

## 1. Coherent-integration data analysis: Rust + Wasm + WebGPU

Open **[juha.no/fmcw/lab/](https://juha.no/fmcw/lab/)**. Select the blue background
interval, select the purple analysis interval, and press **Play**. Results appear
during the scan. Download the completed analysis as **HDF5**.

The shared bare link loads the configured default view: eight chirps, background
[3633,3692), analysis [3704,3880), range 0.001–3 m, radial velocity 0–600 m/s,
and radial acceleration 0–1,000,000 m/s². A `?gui=...` URL overrides defaults
and preserves the controls, selections and zoom.

The matched-filter bank searches midpoint **r0, v0 and a0** using
`R(t) = r0 + v0*(t-t_mid) + 0.5*a0*(t-t_mid)^2`. It applies phase corrections,
then uses fast-time and chirp-time FFTs to accelerate **multi-chirp coherent
integration**. Each receiver has its own background subtraction and noise
estimate; the trajectory search adds the four receiver matched powers
incoherently. The maps show MAX projections over the omitted parameter.
Antenna phase refinement, beamforming, RCS and metallic-sphere diameter estimates
follow the trajectory fit. Scan histories use unconnected points.

**Where the analysis code lives:**

| File | Responsibility |
| --- | --- |
| [`lab-core/src/lib.rs`](lab-core/src/lib.rs) | Rust DSP engine, complex64 FFT, background subtraction, templates and Wasm API |
| [`lab-core/src/radial_search.rs`](lab-core/src/radial_search.rs) | Automatic r0/v0/a0 grid, correction groups, sparse FFT bank, sampling-alias branch retention and joint GLS motion refinement |
| [`web/lab/worker.mjs`](web/lab/worker.mjs) | Loads data, runs searches and incremental analysis scans |
| [`web/lab/gpu-radial.mjs`](web/lab/gpu-radial.mjs) | WebGPU compute shaders for complex64 FFTs and scoring |
| [`web/lab/radial-backend.mjs`](web/lab/radial-backend.mjs) | Backend selection, GPU/CPU consistency checks and CPU fallback |
| [`web/lab/receiver-trains.mjs`](web/lab/receiver-trains.mjs) | Per-RX background residuals, noise powers and receiver covariance |
| [`web/lab/joint-beam.mjs`](web/lab/joint-beam.mjs), [`noise-metric.mjs`](web/lab/noise-metric.mjs) | Analytic complex-amplitude fit and covariance-weighted beamforming; legacy phase-grid routines remain in `beamforming.mjs` for reference |
| [`web/lab/rcs.mjs`](web/lab/rcs.mjs) | Radar equation, PEC-sphere Mie RCS and diameter inversion |
| [`web/lab/app.mjs`](web/lab/app.mjs), [`plots.mjs`](web/lab/plots.mjs) | GUI, graphical windows, WebGL colormesh and waveform/history plots |
| [`web/lab/gui-state.mjs`](web/lab/gui-state.mjs) | Default setup and URL state |
| [`web/lab/export-analysis.mjs`](web/lab/export-analysis.mjs) | HDF5 analysis export |
| [`tools/prepare_lab_data.py`](tools/prepare_lab_data.py) | NetCDF recording → HDF5 and browser float32 transport |

**Run locally**, from the repository root:

```sh
# Rebuild only when changing Rust; core.wasm is already tracked.
rustup target add wasm32-unknown-unknown
./tools/build_lab.sh

# Raw measurement data are not in Git. Import your recording first.
conda run -n base python tools/prepare_lab_data.py /path/to/recording.nc
conda run --no-capture-output -n base python -m http.server 18791 --bind 127.0.0.1 --directory web
```

Open **http://127.0.0.1:18791/lab/**. The DSP runs locally in the browser.
WebGPU accelerates processing where available; Rust/Wasm is the CPU fallback.
WebGL2 renders the colormesh plots. The importer requires NumPy, h5py and netCDF4.

**Tests**, from the repository root:

```sh
cargo test --manifest-path lab-core/Cargo.toml
node tools/validate_complexity_counts.mjs
node tools/validate_gui_state.mjs
node tools/validate_gui_defaults.mjs
node tools/validate_background_residual_noise.mjs
node tools/validate_snr_display.mjs

# Native WebGPU/CPU comparison, fallback and cancellation tests (no browser).
npm install --prefix /tmp/fmcw-webgpu-test webgpu
conda run -n base python tools/run_gpu_validation.py
node tools/validate_gpu_worker.mjs
```

Measurement-dependent tests require prepared Test 63 data. See
[lab instructions](web/lab/README.md), [GPU validation and benchmarks](web/lab/GPU.md),
and [beamforming/RCS conventions](web/lab/BEAMFORMING.md).
The [complexity memo](https://juha.no/fmcw/slides/radial-fft/fmcw_memo_010.pdf)
compares 1, 2, 4, 8, 16 and 32 chirps; regenerate its counts with
`conda run -n base python tools/radial_complexity.py`.

## 2. Measurement simulator

[`web/sim/model.mjs`](web/sim/model.mjs) implements the Equation 25 retarded
FMCW signal model. The simulator shows target motion, chirp timing, complex
beat voltage and a spectrogram. Voltage is relative, not calibrated ADC volts.
It generates simulated signals; it does not analyze recordings.

```sh
conda run --no-capture-output -n base python -m http.server 18765 --bind 127.0.0.1 --directory web
```

Open **http://127.0.0.1:18765/sim/**. No npm build is needed; Three.js is vendored.
Validate against Python with `conda run -n base python web/sim/validate_model.py`.
See [simulator model and units](web/sim/README.md).

## 3. Manim Slides

| Source | Scene | Content |
| --- | --- | --- |
| [`manim/radial_fft_search.py`](manim/radial_fft_search.py) | `RadialFFTSearch` | **42 slides:** plane-wave and complex baseband derivation, coherent integration, animated round-trip propagation phase, full echo phase alignment, FFT matched-filter implementation, receive-array beamforming, SNR → RCS → Mie diameter, and the GUI |
| [`manim/fmcw_space_radar.py`](manim/fmcw_space_radar.py) | `FMCWSpaceRadar` | 47-slide radar lecture, simulations and measured examples |
| [`manim/radar_equation_noise.py`](manim/radar_equation_noise.py) | `RadarEquationNoise` | Radar equation, thermal noise and coherent integration |

The receive-array diagrams and beamforming animations are implemented in
[`manim/radial_beamforming_slides.py`](manim/radial_beamforming_slides.py), using
the PCB layout in [TI’s xWR1843BOOST guide, figure 10](https://www.ti.com/lit/ug/spruim4b/spruim4b.pdf).
[`manim/radial_fft_integration_slides.py`](manim/radial_fft_integration_slides.py)
explains full-vector equivalence, SNR optimality assumptions, and the sparse
two-stage FFT work estimate. Animated grids show why the first FFT computes
all bins while the across-chirp FFT uses the columns requested by the physical
search. Range–Doppler diagrams distinguish coupling from wrapped-frequency
aliases; search bounds alone do not guarantee a unique solution.

Render and export the new deck:

```sh
cd manim
conda run -n base python radial_fft_assets.py
conda run --no-capture-output -n base manim-slides render --disable_caching -r 1920,1080 --fps 30 radial_fft_search.py RadialFFTSearch
conda run -n base python export_web.py qa_fmcw_web RadialFFTSearch
conda run --no-capture-output -n base python -m http.server 18766 --bind 127.0.0.1 --directory qa_fmcw_web
```

Open **http://127.0.0.1:18766/radial-fft/**. Manim needs FFmpeg, LaTeX and dvisvgm;
Python dependencies are in [`requirements.txt`](requirements.txt). The exporter
checks every slide's final decoded frame and installs the final-slide guard.
The new radial deck includes an **Open laboratory** link; detailed derivations
remain in the memos. The older lecture decks retain their provenance footers.
See [slide sources and rendering](manim/README.md).

## Scientific reference code and data

[`capture_analysis/`](capture_analysis/) contains the Python reference signal
model, sphere scattering solver, simulation and offline measurement-analysis
scripts. [`output/`](output/) and [`manim/assets/`](manim/assets/) contain selected
lecture figures, provenance and small HDF5 rendering inputs.

**Keep recordings out of Git.** Deployment data live under `/mnt/shovel/fmcw/`
on `juha-no`; scientific data products use HDF5. `.f32` files are browser/GPU
transport, not archival products. Test 63 timing is reconstructed from editable
assumptions, not measured frame timestamps. RCS/diameter estimates depend on
system-temperature, power, antenna-gain and beamforming assumptions. A fitted
trajectory alone does not establish a detection.

The default radial FFT analyzer uses a joint least-squares final step. It fits
`z_j = q_j(r0,v0,a0) A + n_j`, with one complex amplitude per receiver and
receive covariance estimated from the mean-subtracted quiet raw beat signal.
For every motion trial the amplitudes are solved analytically; bounded
Nelder–Mead adjusts motion using the profiled covariance-weighted residual.
The same fit supplies beam weights proportional to `C^-1 A`. Noise is assumed
white across acquired sample times. The beam signal-power estimate removes the
noise contribution from fitting the complex amplitudes; trajectory search can
still bias weak peaks. RCS uses actual covariance-calibrated combining gain,
assuming equal calibrated physical antenna responses.

Distinct range/velocity candidates retain their identity when FFT indices wrap.
The candidate shortlist includes the strongest representative of every
sampling-alias branch in the configured bounds. Full-waveform checks can reject
approximate aliases, but identical acquired waveforms remain ambiguous; bounds
do not exclude contamination by an echo outside them. Alternative branch checks
are available beside the fit summary and in the HDF5 export.

Mie diameter results retain all inversion roots and also report their min–max
envelope within the chosen diameter bounds. This envelope is not a confidence
interval and does not imply every intermediate diameter is consistent.

Validate the joint fit with `node tools/validate_joint_gls.mjs` and
`node tools/validate_joint_wasm.mjs`; the latter checks actual four-receiver
Rust/Wasm trajectory recovery and the analytic weighted score.

At the default 2 dB full coherent processing-loss setting, the 16-chirp bank contains
28,440,405 nodes. The browser cap is 64 million nodes; candidate selection uses
a bounded heap rather than allocating an index for every node. GPU storage
limits may require CPU fallback. A default 32-chirp bank (~512 million nodes)
still requires narrower bounds or a different storage strategy.

The slide equations and symbol definitions use MathTex or explicit inline TeX math. The least-squares derivation expands the complex residual, identifies the matched inner product, and profiles out the complex amplitude before the trajectory search.
