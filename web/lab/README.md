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
  All three colour plots share one interval: drag the blue Background window to move it or its edges to
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

All three maps label the gold coherent-integration interval. The best-fit
voltage is followed by the continuous geometric slant range
`R(t) = sqrt((v*t-x0)^2 + y0^2)` over the analyzed train. This is inferred
from the fitted trajectory, not an independently measured range.
The complex residual (background-subtracted measurement minus the fitted
complex template, both Re and Im) appears immediately below the voltage fit.

The measurement zoom slider starts at the full recording and zooms around the
coherent train midpoint; selecting another train recenters an active zoom.
Manual bounds and other view buttons retain custom views. Matched-filter maps
always show the complete configured search grid and mark the global best fit
with a cross; colour-scale maxima equal the actual matched-filter peak.
The yellow integration window is draggable on all three measurement maps.
Its pulse count stays fixed and trains stay within one frame. Releasing a
moved train repeats an existing search; otherwise it selects the train.

## Coherent correction-grid FFT search (Memo 8)

The automatic mode is intended for **multi-chirp** integration. It uses
midpoint coordinates `R = range`, `U = radial velocity`, and `G = range
acceleration`. For each velocity branch,
`v = ±sqrt(U² + R*G)`, `x0 = v*t_mid - R*U/v`, and
`y0 = sqrt(G*R³)/abs(v)`. Original human-supplied x0/y0/v bounds remain the
constraints. A separate static template handles v=0. One-chirp comparison
retains the existing direct reference; its ordinary Fourier spectrum is
already available in the measurement plots.

The exact retarded phase is separated into a constant, a fast/slow frequency
plane, and its nonlinear remainder. Similar correction curves share a
zero-padded 2D FFT over fast samples and chirp index. There are no fabricated
samples in chirp gaps and no independent fitted phase per chirp. The fast
frequency mixes range and Doppler; physical slow-frequency aliases are
checked during refinement. Both stored-I/Q orientations remain supported.

FFT energy peaks propose extra physical candidates. The quiet-referenced
statistic, projection points, and fitted waveform are evaluated using the
full physical template at every accepted representative and every refined
candidate. This screening/refinement is approximate and does not guarantee
the global continuous maximum. The maps take MAX over evaluated physical
points binned to the displayed axes; gray map pixels have no evaluation.
They do not fill a physical bounding rectangle with a single point's score.

Grid spacing comes from sampled phase changes/derivatives at feasible local
nodes, with a safety factor. **Estimated cell phase loss is not a certified
worst-case loss.** Correction sharing has a separate tolerance; FFT bins and
finite refinement add approximations. Full bounds close to the radar can
still exceed the 150,000-cell cap, particularly for 8/16 chirps. The worker
reports failure and does not return an incomplete bank or silently shrink
bounds. Narrower physically justified bounds can be much cheaper.

Validate the shipped binary and benchmark 4/8/16-chirp synthetic recovery:

```sh
conda run --no-capture-output -n base python tools/validate_fft_search.py
```

The comparison uses identical measurements and bounds, against a 31×11×41
physical reference grid; these are different discretizations. Timings include
automatic planning, FFT search, and physical refinement, and are specific to
the machine. Exact returned scores, fitted complex waveforms, mirrored geometry,
and reversed stored I/Q are verified independently. Results and the Wasm hash
are stored in `web/lab/qa/fft_search_validation.h5` (ignored by Git).
Implementation: `lab-core/src/fft_search.rs`; benchmark driver:
`tools/benchmark_fft.mjs`; memo: Overleaf project `6ac509d9d83f597ac920e7de`,
`fmcw_memo_008.tex`.
