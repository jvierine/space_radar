# FMCW signal explorer

Public destination: https://juha.no/fmcw/sim/ . Directly accessible, noindex,
and unlisted from the space course index.

Serve this directory over HTTP; ES modules cannot be opened through file://.
For example, from the radar workspace:

```sh
conda run --no-capture-output -n base python -m http.server 18765 --bind 127.0.0.1 --directory web
```

Open http://127.0.0.1:18765/sim/ . All runtime dependencies are vendored;
Three.js 0.179.1 and its MIT license are in vendor/.

## Model and units

`model.mjs` implements the memo Equation 25 phase and retarded reflection
geometry. The radar is at (x0,y0), the target at (v*t,0). x0 is along the
trajectory, not slant range; v is signed along-path velocity, not radial
velocity. Positions are metres, times seconds, frequencies Hz, and velocity
m/s internally. The page supports one to ten chirps and +/-14 km/s.

The reference profiles have 128 samples at 12.5 MS/s (10.24 us live ADC),
0.37 us ramp tail, and 3.10 us idle. Their 13.71 us repetition interval is
inferred from the reader and explicitly labeled an assumption on the page.
The target moves continuously through these gaps. No zero samples are
inserted; curves and STFT windows stay within each acquired chirp.

Complex voltage uses received RF times conjugated transmitted RF. Approaching
carrier Doppler is positive. Amplitude is proportional to 1/R^2 for fixed
cross section, normalized to a return at 1 m. It is not calibrated voltage,
ADC counts, or a hardware feasibility calculation. Analogue frequency and
its ADC alias are displayed separately.

The optional receiver view approximates quasistatic 175/350 kHz high-pass
filters and an ideal -10 to 0 MHz Complex-1x sideband mask. Filter transients
and the true DFE response are not modeled. The ideal view has no filters.
Cross-ramp return samples at the start of a chirp are excluded.

The spectrogram uses 64-sample Hann windows, hop 16, FFT 256. The 195.3 kHz
Fourier scale is set by the 5.12 us window, not the zero-padded bin spacing.
Its dB scale is relative to the peak in the displayed sequence.

Play slows the full sequence to seven seconds and repeats it. The target,
plot cursors, chirp highlight, and ADC/tail/idle labels share the same time.
Scrubbing pauses playback. The 3D animation shows the moving target against
the static radar and trajectory geometry. Light-time delay is computed in the
signal model.

`chip.svg` is a simplified signal-flow diagram based on the TI IWR1843
datasheet functional block diagram, not a wiring or pin diagram.

## Verification

```sh
conda run -n base python web/sim/validate_model.py
```

This compares 16 profile/geometry/chirp-count cases with the existing memo
Python implementation, checks padding and global sample time, differentiates
phase against analytic IF, and checks FFT signs and aliasing with known tones.
Results and source hashes are saved to `qa/model_validation.h5`.

Browser verification covers both profiles, both receiver modes, +/-14 km/s,
zero velocity, one and ten chirps, geometry fields, play/pause, and scrubbed
ADC/ramp-tail/idle/sequence-end states. Visual verification includes rendered
plots after initial layout resizing. No package or environment changes are
required.

Deployment copies only index.html, style.css, app.mjs, model.mjs, chip.svg
and vendor/ to /var/www/html/fmcw/sim/ on juha-no. QA files and captures are
local products, not public website assets.
