# Matched filtering and incremental analysis

The trajectory bank sums **powers from all four receivers**, after subtracting
each receiver's own complex mean at every fast-time sample. At each template
and common stored-I/Q orientation the score is

```
rho = sum_rx |MF(event_rx)|^2 / mean_quiet_train(sum_rx |MF(quiet_rx)|^2)
display = 10 log10(rho)
```

Quiet trains are common to the four receivers, nonoverlapping, frame-contained,
and outside both event and background-mean intervals. Padding/nonfinite chirps
on any receiver exclude that train. The display receiver selector affects raw
plots and fit traces, not the bank's objective. Beamforming subsequently searches
receiver phases at the shared trajectory fit.

CPU and WebGPU FFT arithmetic/storage are complex64 (two float32 components).
CPU twiddles are cached; phase construction and CPU energy accumulation remain
float64. WebGPU uses WGSL workgroup-memory fast-time and slow-time FFTs, sparse frequency
columns, and the same correction phases/bin plans exported by Rust. Phase plans
and the exact quadratic-template refinement use float64 Rust calculations.
MAX projections remain in the frontend. A completed GPU bank is checked against
CPU groups at the beginning, middle and end before import. Adapter, shader,
allocation, device-loss, timeout, or numerical failures discard the bank and
restart it on CPU. Cancellation does not start a fallback search. Small banks
use CPU automatically; CPU can also be selected explicitly.

The purple analysis window specifies first and last analyzed chirps (the stop
control is exclusive). Every contained intact train start is searched, except
quiet intervals and frame crossings. Each completed point immediately updates
the trajectory maps, fit, phase maps and line plots. Receiver-pair angles are
**phase differences**, wrapped to −180…180 degrees; no antenna geometry is
assumed to convert these into arrival directions. Scan histories retain scalar
results, not full cubes. Cancellation retains completed points.

GUI controls and heatmap viewports are encoded in the `gui` URL parameter.
Opening a shared fitted view restores the controls and recomputes its selected
train. Recordings and score arrays are not embedded in the URL. MathJax 3.2.2 SVG
axis labels are self-hosted and included in PNG exports. The circular UiT seal
is the official asset from
https://uit.no/ressurs/uit/profil2019/logo/UiT_Segl_Bok_Bl%C3%A5_RGB.svg.

## Validation

Native Dawn/Metal tests run the shipped WebGPU shaders without a browser:

```sh
npm install --prefix /tmp/fmcw-webgpu-test webgpu
conda run -n base python tools/run_gpu_validation.py
node tools/validate_gpu_worker.mjs
node tools/validate_receiver_sum.mjs
node tools/validate_scan_worker.mjs
node tools/validate_gui_state.mjs
node tools/validate_beam_autorun.mjs
PATH=/Users/jvi019/.cargo/bin:$PATH cargo test --manifest-path lab-core/Cargo.toml
```

The HDF5 benchmark report is `web/lab/qa/gpu_validation.h5` (ignored by Git).
It records source hashes, grids, score errors and cold/warm timings. Full cubes,
valid masks, both MAX projections and refined peaks are compared with CPU.
Tests include conjugated I/Q, receiver filtering, 4/8/16 chirps, opposing receiver
phases, unequal receiver noise, independent backgrounds, cancellation with
immediate replacement, device loss and injected fallback failures.

Test 63, train 750, four receivers, eight chirps and the complete default
18,551,652-node bank: **43.77 s complex64 CPU versus 4.10 s warm WebGPU (10.67×)**;
cold WebGPU 3.99 s. Timings include CPU consistency checks and exact peak
refinement. Maximum normalized score error was 1.58e−6. These are native
Dawn/Metal measurements on this Mac, not browser timings. An actual browser
can have different dispatch overhead and hardware/limit availability.
