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
those identical reference blocks. Three phase maps show RX1 × RX2 (MAX over
RX3), RX1 × RX3 (MAX over RX2), and RX2 × RX3 (MAX over RX1), in one row.
All use the same color limits, with the actual grid peak as their maximum,
and crosses at the projected best grid cell. Nelder–Mead refines the reported
phases and peak; the maps show the discrete phase grid.

The receiver-growth table adds RX1, RX2, and RX3 to RX0 in that order,
optimizing the phases separately for each subset. All use the same quiet
trains. Ideal added SNR equals the added receiver's individual linear ratio
under independent noise and optimal amplitude weights. Percentage is
`100 * (subset_ratio - previous_subset_ratio) / added_receiver_ratio`.
Phase-only equal-amplitude combining can reduce the result; correlated-noise
cancellation can exceed this independent-noise benchmark.

Total gain uses the arithmetic mean of the four individual matched-energy
ratios in **linear** units as its single-antenna baseline. Divide the refined
four-receiver ratio by that baseline; ideal total gain is 4× (6.02 dB), and
percentage of ideal total SNR is `100 * total_gain / 4`. This total baseline
differs from the ordered receiver-growth table's RX0 baseline. The total
percentage compares full linear gain; the per-addition percentage compares
the added linear SNR relative to the added receiver alone.

## Conditional RCS and sphere diameter

Use 9000 K by default, with editable 12 dBm TX power, 6 dBi TX/RX gains,
extra loss, and distance initialized from the fitted midpoint range. These
power/gain defaults follow the existing capture-analysis assumptions.
At the fast-time midpoint frequency, infer received power from
`max(matched_ratio - 1, 0) * k*T / effective_time`, where
`effective_time = (sum(abs(q)))^2 / (fs * sum(abs(q)^2))`.
The monostatic radar equation gives sigma, using constant range/RCS over
the train. The beamformed row divides by ideal four-RX receive gain.
Actual correlated quiet noise and selection bias are not calibrated away.

`rcs.mjs` uses the PEC monostatic Mie series from the existing
`capture_analysis/analyze_captures.py:sphere` convention. Inversion lists
all sign-crossing roots on a 4096-interval logarithmic diameter grid within
the editable bounds, followed by bisection. Multiple solutions are retained;
an exact tangent root can require tighter diameter bounds. Run
`conda run -n base python tools/validate_rcs.py` for comparison to SciPy,
multiple-root inversion, noise-floor subtraction, and beam gain scaling.

The reported gain is in-sample, after optimizing phase on the event and quiet
reference. Searching many phases raises noise maxima. It is not independent
proof of detection or a calibrated thermal SNR. Antenna positions, gain/phase
calibration and a physical steering model are needed to infer an arrival angle.

Run `node tools/validate_beamforming.mjs`. Checks cover known phase recovery,
the 6.02 dB gain for four equal signals with independent receiver noise,
fully correlated noise without fictitious array gain, separate complex means,
and packet padding on a single receiver.
