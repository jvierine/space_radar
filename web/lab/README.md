# FMCW laboratory

Live browser: https://juha.no/fmcw/lab/ . Unlisted from the space course index.
Data on the deployment host: `/mnt/shovel/fmcw/`.

## Build and run

The Rust engine has no external crate dependencies and exposes a small C ABI
to JavaScript. It runs inside a Web Worker; WebGL2 renders float textures.
The DSP runs on the user's browser, with no remote analysis service.

```sh
rustup target add wasm32-unknown-unknown
./tools/build_lab.sh
cargo test --manifest-path lab-core/Cargo.toml
conda run --no-capture-output -n base python tools/prepare_lab_data.py /path/to/recording.nc
conda run --no-capture-output -n base python tools/validate_lab.py
conda run --no-capture-output -n base python -m http.server 18791 --bind 127.0.0.1 --directory web
```

Open http://127.0.0.1:18791/lab/ . Python needs NumPy, h5py and netCDF4 for
import and validation. Node.js runs the compiled-Wasm validation harness.
The importer writes `recording.h5` as the scientific product and little-endian
interleaved float32 `.f32` arrays as lossless HTTP/GPU transport. It verifies
every transported complex sample against the original. The browser checks
the selected RX transport SHA-256 before importing it into Wasm memory.
Prepared data and QA products are ignored by Git; raw captures are not bundled.

## Plots and timing

- Raw means the stored complex ADC signal after analogue chirp downconversion,
  with no FFT or matched decoding. Select its real or imaginary component.
- The quiet interval is user-selected, not an automatically identified trigger.
  A complex mean is formed per fast-time sample over intact quiet chirps. The
  Spectrum plot: drag the blue Background window to move it or its edges to
  resize; release to recompute. Focus its buttons and use arrow keys (Shift
  for ten chirps) for precise adjustment. The
  SAME mean is subtracted from all displayed and analyzed chirps.
- FFTs are complex, Hann-windowed and zero-padded from 225 to 256 points for
  Test 63. Actual Fourier resolution is fs/225 = 55.56 kHz, while padded bins
  are 48.83 kHz apart. Frequency is signed stored-I/Q beat frequency, which
  combines delay and Doppler. Original and subtracted views share one raw
  recording spectral peak reference.
- The source clock is 2026-07-23 11:02:58; its timezone is unspecified. There
  are no per-frame timestamps. ADC sampling starts 2 µs into the ramp; 225
  samples at 12.5 MS/s collect 18 µs. Idle is 5 µs. Ramp tail is not stored:
  default 0.37 µs is an assumption, giving a 25.37 µs chirp period. Default
  frame interval assumes continuous frames, 125 × 25.37 µs = 3.17125 ms.
  Both are editable. Clock labels are Unix epoch seconds with an explicit
  UTC assumption because the source timezone is unspecified. View bounds use
  inclusive first/last chirps; Zoom to background shows the selected mean interval.
  Clock labels are reconstructed, and concatenated-index
  plots omit any configured gaps. Coherent searches do not cross frames.
- Packet padding is recorded globally without a per-sample mask. A run of
  at least eight exact complex zeros, or a nonfinite sample, flags a chirp.
  Flagged columns are grey; traces retain original samples for inspection.
  Whole flagged chirps are excluded from mean, spectra and filter statistics.
  Isolated zero codes are retained. This is a conservative QC inference.

## Matched filter

The exact retarded Eq. 25 phase is
`2π[-f0 τ − S u τ + S τ²/2]`. Global time within a train uses actual sampled
positions with ADC/idle gaps. The echo voltage amplitude is proportional to
`1/R_ret²`, assuming constant scattering strength. Initial radar coordinates
are `(x0,y0)`, the projectile position is `(v*t,0)`, and `x0` is longitudinal
offset, not initial slant range. Each selected train is phase coherent.
Optional quasistatic 175/350 kHz high-pass factors approximate receiver
filtering; their configuration and transient response are not known for this
recording. Neither a beam gate nor a particular one-sided I/Q convention is
silently imposed. Both I/Q orientations are searched.

For each template q and residual data z, the score is
`rho = |qᴴ z|² / mean_quiet(|qᴴ n|²)`, reported as `10 log10(rho)`.
The template norm cancels. The denominator is estimated separately for that
template and I/Q orientation from up to 24 evenly selected, nonoverlapping,
intact quiet-reference trains of the same length. Those trains cannot overlap
the selected train or cross a frame boundary; at least eight are required.
The default mean and noise-reference ranges are disjoint. If the selected
train overlaps the mean interval, the interface identifies that overlap.

This measures observed matched output relative to residual quiet-background
noise/clutter, not isolated signal energy or calibrated thermal SNR. It is
not an optimally covariance-whitened detector. A single fixed noise-only
filter has mean linear rho=1 under a stationary noise model; searching many
correlated filters raises the maximum. There is no universal 13 dB threshold.
The short-train x0/y0/v ambiguities remain visible in the maximum projections;
mirroring `(x0,v)` to `(-x0,-v)` preserves monostatic range history.

The time scan tests valid starts in a configurable interval and stride,
excluding the mean and noise-reference intervals. A scan is limited to 512
starts and a bank to 150000 templates; use a coarse grid/stride and then refine.
After a scan, maps and the fitted trace refer to its highest-scoring start.
Comparison searches show all five pulse counts and display the best score's
maps; changing numerical settings invalidates previous maps.

## Validation and provenance

`tools/validate_lab.py` exercises the compiled browser binary, not a JavaScript
replacement. It compares 40 retarded templates with the memo model, validates
the quiet complex mean and complex Hann FFT against NumPy, compares complete
1/2/4/8/16-pulse score cubes, checks an injected coherent waveform/time, and
rejects padding and frame-crossing trains. Synthetic injection does not prove
that the measured recording contains a detected projectile or that short
trains uniquely recover its trajectory. Results/hashes are saved in HDF5.

Rust source: `lab-core/src/lib.rs`; GPU/axis plots: `web/lab/plots.mjs`;
worker orchestration: `web/lab/worker.mjs`; controls: `web/lab/app.mjs`;
original data import: `tools/prepare_lab_data.py`.

Platform references: [Rust wasm32 target](https://doc.rust-lang.org/rustc/platform-support/wasm32-unknown-unknown.html),
[WebGL texture upload](https://developer.mozilla.org/en-US/docs/Web/API/WebGLRenderingContext/texImage2D).
