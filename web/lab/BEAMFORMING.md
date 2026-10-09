# Four-receiver phase search

Fit a trajectory, then use **Search 1,000 phase combinations** below the fit.
RX0 fixes the phase gauge at zero. RX1, RX2 and RX3 each use ten phases
0, 36, ..., 324 degrees; phase steps are adjustable from 2 to 32.
The combined voltage is `(RX0 + exp(iφ1) RX1 + exp(iφ2) RX2 + exp(iφ3) RX3)/2`.
It is phase-only beamforming with fixed equal amplitudes, conditional on the
selected trajectory and stored-I/Q orientation. It does not change the
single-receiver measurement panels or perform a joint trajectory search.

Each stream is length- and SHA-256-verified. The background mean is estimated
separately for each receiver using chirps valid on all four receivers. Padding
or a frame-crossing selected train is rejected. Up to 24 common independent
quiet trains outside both the selected train and mean-estimation interval
supply the reference; at least eight are required.

The Rust/Wasm Eq.25 template supplies per-channel complex matched projections.
Every phase combination is applied to the event projections and all quiet
projections. Score is event matched energy divided by average quiet matched
energy, retaining cross-receiver noise correlation. Single-channel scores use
those identical reference blocks. The phase plot is MAX over RX3, with the
actual peak as its color-scale maximum and a cross at the best cell.

The reported gain is in-sample, after optimizing phase on the event and quiet
reference. Searching many phases raises noise maxima. It is not independent
proof of detection or a calibrated thermal SNR. Antenna positions, gain/phase
calibration and a physical steering model are needed to infer an arrival angle.

Run `node tools/validate_beamforming.mjs`. Checks cover known phase recovery,
the 6.02 dB gain for four equal signals with independent receiver noise,
fully correlated noise without fictitious array gain, separate complex means,
and packet padding on a single receiver.
