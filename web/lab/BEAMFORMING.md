# Joint motion and receiver-amplitude fit

The default radial FFT search first uses incoherently added receiver powers to
locate candidate trajectories. The final bounded Nelder–Mead fit minimizes
`sum_j (z_j-q_j A)^H C^-1 (z_j-q_j A)`, where `z_j` contains four complex RX
samples, `q_j` is the complete quadratic-range FMCW echo, and `C` is the
quiet raw residual receive covariance. Temporal noise is assumed white.
The four complex amplitudes are solved analytically as `A=s/Q`, with
`s=sum_j conj(q_j) z_j` and `Q=sum_j |q_j|^2`. Thus the optimizer maximizes
`s^H C^-1 s/Q`. The same fit supplies complex beam weights proportional to
`C^-1 A`; there is no separate phase-only grid or phase simplex in this path.

The covariance inverse uses Cholesky factorization after normalizing the RX
variances. A small diagonal correlation ridge is used only if needed for
positive definiteness, and is reported/exported. Beam noise is propagated
through the actual complex weights, including cross-RX correlations.

Fitting four complex amplitudes also fits noise. The GLS quadratic score has
noise contribution `tr(C_reg^-1 C)`, normally four. Subtract this contribution
before converting fitted signal SNR to received power. Further trajectory
search selection bias remains. Displayed dB values have the requested 0 dB
floor; physical estimates use unfloored linear signal power.

RCS combining gain uses the implemented weights and measured receive
covariance, assuming equal calibrated antenna responses in each RX's noise
units. It equals four for equal independent RXs, but is not assumed universally
four. Fitted receiver phases include electronic offsets and are not calibrated
arrival angles. HDF5 stores complex amplitudes, weights, covariance, GLS score,
noise bias and calibration gain.

## Reference implementation retained below

The old phase-only routines in `beamforming.mjs` remain as tested reference
utilities, not the active default radial processing path.

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

After a trajectory search completes, automatically process all four receivers,
refine the beam phases, and display the RCS/diameter table above the three
phase marginals. The phase-search button can rerun this calculation.

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
