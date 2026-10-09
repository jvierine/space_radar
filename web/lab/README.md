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
  Both are editable. Time labels show reconstructed seconds since file start,
  independent of the source clock timezone. View bounds use
  inclusive first/last chirps; Zoom to background shows the selected mean interval.
  Clock labels are reconstructed, and concatenated-index
  plots omit any configured gaps. Coherent searches do not cross frames.
- Packet padding is recorded globally without a per-sample mask. A run of
  at least eight exact complex zeros, or a nonfinite sample, flags a chirp.
  Flagged columns are grey; traces retain original samples for inspection.
  Whole flagged chirps are excluded from mean, spectra and filter statistics.
  Isolated zero codes are retained. This is a conservative QC inference.

## Matched filter

The direct geometry reference uses the exact retarded Eq. 25 phase
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
starts and the direct geometry bank to 150000 templates; use a coarse grid/stride and then refine.
After a scan, maps and the fitted trace refer to its highest-scoring start.
Direct comparison searches show all five pulse counts and display the best score's
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
voltage is followed by the inferred range over the analyzed train: the
quadratic polynomial in radial mode, or `sqrt((v*t-x0)^2+y0^2)` in the direct
geometry reference. Neither is an independently measured range.
The complex residual (background-subtracted measurement minus the fitted
complex template, both Re and Im) appears immediately below the voltage fit.

The measurement zoom slider starts at the full recording and zooms around the
coherent train midpoint; selecting another train recenters an active zoom.
Drag a rectangle in empty colormesh space to zoom. Measurement plots share
the horizontal chirp interval; voltage plots share the fast-time interval.
Background bodies, labels, and resize edges and coherent-window bodies and
labels keep their existing drag actions. Reset full view restores the recording
and vertical axes. Matched-filter and antenna-phase maps have local reset buttons.

The antenna-phase search refines its best discrete cell with Nelder–Mead over
the three relative receiver phases, with RX0 fixed at zero. Both stages use
the same background-referenced matched-energy ratio and quiet trains. The
reported phases and peak are refined; the projection remains the discrete
grid, and the table reports both peaks.

The viewer uses a white background, dark plot labels, and brown/teal traces.
Each plot's PNG button exports its current view at four times the displayed
pixel dimensions, with freshly rendered labels, axes, and color bar on an
opaque white background. Interactive controls are excluded from the figure.

Beamforming includes three MAX phase projections, separate subset phase
optimization for receiver-growth comparisons, and total gain relative to the
linear mean of the four individual ratios. Expand the receiver diagnostics
to compare each channel's observed and quiet matched powers. Background
means are removed separately for each channel using common intact chirps.
The RCS/PEC-diameter table uses editable thermal/radar assumptions (9000 K
by default); see `BEAMFORMING.md` for its conversion and inversion conventions.

Manual bounds and other view buttons retain custom views. Matched-filter maps
initially show the complete configured search grid and mark its global peak
with a cross; radial FFT screening peaks and directly verified fits are
labeled separately. Colour maxima equal the actual plotted grid peak.
The yellow integration window is draggable on all three measurement maps.
Its pulse count stays fixed and trains stay within one frame. Releasing a
moved train repeats an existing search; otherwise it selects the train.

## Explicit radial FFT search (Memo 9)

Automatic mode searches independent **midpoint r0, v0, a0 bounds**. The model
is `r(t)=r0+v0*h+0.5*a0*h²`, with `h=t-t_mid`. r0 is one-way delay range;
v0 is signed radial velocity (positive receding); a0 is radial acceleration
in m/s². Its carrier Doppler chirp rate is approximately `-2*f0*a0/c` Hz/s.
The polynomial is a local approximation, with cubic and higher range terms
omitted. It is not the geometry reference's Cartesian speed/acceleration.
The radial template is phase-only with a complex fitted gain; the direct
geometry template retains its `1/R_ret²` envelope. Optional HPFs affect
radial candidate verification, while FFT screening is phase-only.

As in [hard_target](https://github.com/jvierine/hard_target/blob/main/src/hardtarget/analysis/gmf/gmf_cpu_numpy.py),
loop through acceleration corrections, then FFT to search velocity. FMCW adds
range beat and a mixed fast/slow-time correction that depends on velocity.
The implementation shares corrections across bounded velocity groups and
computes a sparse 2D FFT: each chirp's fast FFT, then only the requested fast
frequency columns' slow FFTs. Actual measured sample clocks are used without
inventing observations during idle gaps. Both I/Q orientations and all slow
frequency aliases within the supplied velocity bounds are searched.

Acceleration grid spacing and velocity-sharing widths use analytic global
bounds on derivatives of the polynomial phase remainder. The requested loss
allowance sets `epsilon=acos(sqrt(1-loss))` separately for nearest-acceleration
mismatch and correction sharing. These are separate stagewise allowances;
FFT-bin rounding, range/velocity sampling, and polynomial truncation add
errors. It is not a guarantee of total energy loss or global recovery.
Range and velocity axes have padded FFT-derived spacing; all axes include
the original bounds and support singleton intervals. Signed a0 is allowed.
Eight-chirp default bounds r0=[0.001,3] m, v0=[-900,900] m/s and a0=[0,1e6]
m/s² produce 329×1524×37 radial nodes with 555 shared correction FFTs at 5%.
A 32-million-node cap reports an error instead of shrinking the bounds.

Every grid node with positive range throughout the train gets an approximate
FFT matched-energy/quiet-energy ratio. The same correction and bins are
applied to each of the up to 24 independent quiet-reference trains. Strongest
128 candidates are checked by direct inner products with the actual quadratic
template in both orientations. The reported fitted voltage and score are
directly verified and followed by bounded Nelder–Mead refinement in the
continuous r0/v0/a0 domain; this local refinement does not guarantee its global
maximum. Maps are complete MAX projections v0×r0 (over a0) and v0×a0 (over r0).
Crosses and colour maxima refer to the actual FFT grid peak, which is labeled
separately from the verified fit and may be a different node.

The radial mode defaults to eight chirps and compares 2/4/8/16 chirps. The
one-chirp spectrum is already provided above; one-chirp matched searches
remain available only in the existing direct geometry reference. Switching
methods preserves each method's own bounds. Radial time scans apply the same
bank independently to every valid start and show the best verified candidate;
every contained intact start is analyzed. The purple interval controls its
bounds; maps, phase estimates and scalar history plots update point by point.

Validate with:

```sh
conda run --no-capture-output -n base python tools/validate_radial_search.py
```

This exercises the compiled Wasm default 8-chirp domain, reversed I/Q with
aliased radial velocities, and focused 16-chirp searches. NumPy independently
checks the full quadratic phase, selected quiet controls, exact best-score
ratio and fitted voltage. Synthetic data, projections and the binary hash are
saved in `web/lab/qa/radial_search_validation.h5`, ignored by Git. It does not
establish a measured detection or unique trajectory. Timings are specific to
the machine and browser. Source: `lab-core/src/radial_search.rs`; harness:
`tools/benchmark_radial.mjs`; memo: `fmcw_memo_009.tex` in Overleaf project
`6ac509d9d83f597ac920e7de`. The older geometry-coordinate FFT implementation
and its regression harness remain in source but are not the viewer's radial
search path.

See [GPU.md](GPU.md) for four-receiver incoherent scoring, complex64 CPU/WebGPU
FFTs, fallback behavior, URL state, scan outputs and validation commands.
