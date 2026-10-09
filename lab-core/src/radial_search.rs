//! Explicit midpoint polynomial range/range-rate/acceleration bank.
//! Acceleration phasors followed by FFT velocity search, as in hard_target.
use super::*;

pub(super) struct Radial {
    pub axes: [Vec<f64>; 3], // r0, a0, v0 (cube layout)
    pub nf: usize,
    pub ns: usize,
    pub groups: Vec<(usize, usize, usize)>, // acceleration index, velocity begin/end
    pub epsilon: f64,
    pub bin_loss: f64,
    pub refined: Option<[f64; 3]>,
}
fn axis(lo: f64, hi: f64, step: f64) -> Vec<f64> {
    let intervals = ((hi - lo) / step).ceil().max(0.) as usize;
    (0..=intervals)
        .map(|i| {
            if intervals == 0 {
                lo
            } else {
                lo + (hi - lo) * i as f64 / intervals as f64
            }
        })
        .collect()
}
pub(super) fn centre(a: &Engine) -> f64 {
    a.adc + (a.n - 1) as f64 / (2. * a.fs) + (a.pulses - 1) as f64 * a.period / 2.
}
pub(super) fn frequencies(a: &Engine, th: [f64; 3]) -> [f64; 2] {
    let [r, v, _] = th;
    let mid = a.adc + (a.n - 1) as f64 / (2. * a.fs);
    let slow = (-2. * (a.f0 + a.slope * mid) / C + 4. * a.slope * r / (C * C)) * v;
    [slow - 2. * a.slope * r / C, slow]
}
pub(super) fn values(a: &Engine, th: [f64; 3]) -> Vec<Z> {
    let [r0, v0, a0] = th;
    let star = centre(a);
    (0..a.n * a.pulses)
        .map(|i| {
            let u = a.adc + (i % a.n) as f64 / a.fs;
            let h = (i / a.n) as f64 * a.period + u - star;
            let r = r0 + v0 * h + 0.5 * a0 * h * h;
            let tau = 2. * r / C;
            if r <= 0. || u < tau {
                return Z::default();
            }
            let phase = 2. * PI * (-(a.f0 + a.slope * u) * tau + 0.5 * a.slope * tau * tau);
            let mut q = Z {
                r: phase.cos(),
                i: phase.sin(),
            };
            if a.receiver {
                let f = -a.slope * tau - (a.f0 + a.slope * (u - tau)) * 2. * (v0 + a0 * h) / C;
                q = q.mul(Z::hp(f, 175e3)).mul(Z::hp(f, 350e3));
            }
            q
        })
        .collect()
}
fn valid(a: &Engine, th: [f64; 3]) -> bool {
    let [r, v, g] = th;
    let h = (a.pulses - 1) as f64 * a.period / 2. + (a.n - 1) as f64 / (2. * a.fs);
    let mut min = (r - v * h + 0.5 * g * h * h).min(r + v * h + 0.5 * g * h * h);
    if g > 0. && (-v / g).abs() < h {
        min = min.min(r - v * v / (2. * g));
    }
    min > 0.
}
impl Radial {
    pub fn new(a: &Engine, b: [[f64; 2]; 3], loss: f64, cap: usize) -> Option<Self> {
        if b.iter()
            .any(|b| !b[0].is_finite() || !b[1].is_finite() || b[0] > b[1])
            || b[0][0] <= 0.
            || !(loss > 0. && loss <= 0.5)
            || a.slope == 0.
        {
            return None;
        }
        let epsilon = (1. - loss).sqrt().acos();
        let nf = (8 * a.n).next_power_of_two();
        let ns = if a.pulses == 1 {
            1
        } else {
            (8 * a.pulses).next_power_of_two()
        };
        let h = (a.pulses - 1) as f64 * a.period / 2. + (a.n - 1) as f64 / (2. * a.fs);
        let eta = (a.n - 1) as f64 / (2. * a.fs);
        let f = a.f0 + a.slope * (a.adc + eta);
        let k = 4. * PI / C;
        let l = 4. * PI * a.slope / (C * C);
        let vmax = b[1][0].abs().max(b[1][1].abs());
        let amax = b[2][0].abs().max(b[2][1].abs());
        // Global polynomial derivative bound, including FMCW cross terms and tau².
        let da = 0.5 * (k * f.abs() + 2. * l.abs() * b[0][1]) * h * h
            + 0.5 * k * a.slope.abs() * eta * h * h
            + l.abs() * (vmax * h.powi(3) + 0.5 * amax * h.powi(4));
        let dv = k * a.slope.abs() * eta * h + l.abs() * (2. * vmax * h * h + amax * h.powi(3));
        let astep = 2. * epsilon / da.max(1e-30);
        let dr = l.abs() * amax * h * h * (b[0][1] - b[0][0]) / 2.;
        if dr >= epsilon {
            return None;
        }
        let vgroup = 2. * (epsilon - dr) / dv.max(1e-30);
        let rstep = C * a.fs / (2. * a.slope.abs() * nf as f64);
        let vstep = if ns == 1 {
            C * a.fs / (2. * f.abs() * nf as f64)
        } else {
            C / (2. * f.abs() * ns as f64 * a.period)
        };
        // Validate counts before allocation; never shrink the user's bounds.
        let count = |lo: f64, hi: f64, step: f64| -> Option<usize> {
            if !step.is_finite() || step <= 0. {
                return None;
            }
            let intervals = ((hi - lo) / step).ceil();
            if !intervals.is_finite() || intervals < 0. || intervals >= cap as f64 {
                return None;
            }
            (intervals as usize).checked_add(1)
        };
        let nr = count(b[0][0], b[0][1], rstep)?;
        let nv = count(b[1][0], b[1][1], vstep)?;
        let na = count(b[2][0], b[2][1], astep)?;
        if nr.checked_mul(nv)?.checked_mul(na)? > cap {
            return None;
        }
        let axes = [
            axis(b[0][0], b[0][1], rstep),
            axis(b[2][0], b[2][1], astep),
            axis(b[1][0], b[1][1], vstep),
        ];
        let mut groups = vec![];
        for ia in 0..axes[1].len() {
            let mut first = 0;
            while first < axes[2].len() {
                let mut last = first + 1;
                while last < axes[2].len() && axes[2][last] - axes[2][first] <= vgroup {
                    last += 1;
                }
                groups.push((ia, first, last));
                first = last;
            }
        }
        let sinc = |x: f64| if x == 0. { 1. } else { x.sin() / x };
        let bin_loss = 1.
            - sinc(PI * a.n as f64 / (2. * nf as f64)).powi(2)
                * if ns == 1 {
                    1.
                } else {
                    sinc(PI * a.pulses as f64 / (2. * ns as f64)).powi(2)
                };
        Some(Self {
            axes,
            nf,
            ns,
            groups,
            epsilon,
            bin_loss,
            refined: None,
        })
    }
    pub fn theta(&self, i: usize) -> [f64; 3] {
        let na = self.axes[1].len();
        let nv = self.axes[2].len();
        [
            self.axes[0][i / (na * nv)],
            self.axes[2][i % nv],
            self.axes[1][(i / nv) % na],
        ]
    }
    pub fn len(&self) -> usize {
        self.axes.iter().map(Vec::len).product()
    }
}
fn bin(f: f64, rate: f64, n: usize) -> usize {
    ((f / rate * n as f64).round() as i64).rem_euclid(n as i64) as usize
}
// Only transform fast-frequency columns actually requested by the r0/v0 grid.
fn spectrum(
    a: &Engine,
    z: &[Z],
    corr: &[Z],
    conj: bool,
    nf: usize,
    ns: usize,
    columns: &[usize],
) -> Vec<Z32> {
    let mut rows = vec![Z32::default(); a.pulses * nf];
    for k in 0..a.pulses {
        for j in 0..a.n {
            let q = corr[k * a.n + j];
            rows[k * nf + j] = Z32::from(z[k * a.n + j]).mul(Z32 {
                r: q.r as f32,
                i: (if conj { q.i } else { -q.i }) as f32,
            });
        }
        fft32(&mut rows[k * nf..(k + 1) * nf]);
    }
    let mut out = vec![Z32::default(); columns.len() * ns];
    for (i, &j) in columns.iter().enumerate() {
        let col = &mut out[i * ns..(i + 1) * ns];
        for k in 0..a.pulses {
            col[k] = rows[k * nf + j];
        }
        if ns > 1 {
            fft32(col);
        }
    }
    out
}
fn group_plan(
    a: &Engine,
    index: usize,
    orientation: usize,
) -> (Vec<Z>, Vec<usize>, Vec<(usize, usize)>) {
    let bank = a.radial.as_ref().unwrap();
    let (ia, first, last) = bank.groups[index];
    let (nf, ns) = (bank.nf, bank.ns);
    let r = bank.axes[0][bank.axes[0].len() / 2];
    let v = (bank.axes[2][first] + bank.axes[2][last - 1]) / 2.;
    let th = [r, v, bank.axes[1][ia]];
    let mid = a.adc + (a.n - 1) as f64 / (2. * a.fs);
    let f = frequencies(a, th);
    // Phase-only correction stays defined even where the representative range
    // crosses zero; invalid physical grid nodes themselves are masked below.
    let corr: Vec<Z> = (0..a.n * a.pulses)
        .map(|i| {
            let u = a.adc + (i % a.n) as f64 / a.fs;
            let d = ((i / a.n) as f64 - (a.pulses - 1) as f64 / 2.) * a.period;
            let h = d + u - mid;
            let rr = r + v * h + 0.5 * th[2] * h * h;
            let phase = -4. * PI * (a.f0 + a.slope * u) * rr / C
                + 4. * PI * a.slope * rr * rr / (C * C)
                - 2. * PI * (f[0] * (u - mid) + f[1] * d);
            Z {
                r: phase.cos(),
                i: phase.sin(),
            }
        })
        .collect();
    let nr = bank.axes[0].len();
    let na = bank.axes[1].len();
    let nv = bank.axes[2].len();
    let mut nodes = vec![];
    for ir in 0..nr {
        for iv in first..last {
            let i = (ir * na + ia) * nv + iv;
            let th = bank.theta(i);
            if valid(a, th) {
                nodes.push((i, frequencies(a, th)));
            }
        }
    }
    let sign = if orientation == 0 { 1. } else { -1. };
    let mut columns: Vec<usize> = nodes
        .iter()
        .map(|(_, f)| bin(sign * f[0], a.fs, nf))
        .collect();
    columns.sort_unstable();
    columns.dedup();
    let lookup: Vec<(usize, usize)> = nodes
        .iter()
        .map(|(i, f)| {
            (
                *i,
                columns.binary_search(&bin(sign * f[0], a.fs, nf)).unwrap() * ns
                    + bin(sign * f[1], 1. / a.period, ns),
            )
        })
        .collect();
    (corr, columns, lookup)
}

pub(super) fn search_group(a: &mut Engine, index: usize) {
    let (nf, ns) = {
        let b = a.radial.as_ref().unwrap();
        (b.nf, b.ns)
    };
    for orientation in 0..2 {
        let (corr, columns, lookup) = group_plan(a, index, orientation);
        let events: Vec<&Vec<Z>> = if a.events.is_empty() {
            vec![&a.event]
        } else {
            a.events.iter().collect()
        };
        let mut event = vec![0.; columns.len() * ns];
        for z in events {
            let spec = spectrum(a, z, &corr, orientation == 1, nf, ns, &columns);
            for (d, z) in event.iter_mut().zip(spec) {
                *d += z.power();
            }
        }
        let mut noise = vec![0.; event.len()];
        for z in &a.noise {
            let spec = spectrum(a, z, &corr, orientation == 1, nf, ns, &columns);
            for (d, z) in noise.iter_mut().zip(spec) {
                *d += z.power() / (a.noise.len() / a.channels()) as f64;
            }
        }
        for (i, j) in lookup {
            if noise[j] > 1e-24 {
                let score = (event[j] / noise[j]) as f32;
                if !a.scores[i].is_finite() || score > a.scores[i] {
                    a.scores[i] = score;
                    a.signs[i] = orientation as f32;
                }
            }
        }
    }
}
// Float32 complex arrays; index fields are u32 bit patterns, never rounded floats.
// Header: samples, pulses, nf, ns, quiet trains, nodes, columns.
pub(super) fn gpu_group(a: &mut Engine, index: usize, orientation: usize) {
    let (corr, columns, lookup) = group_plan(a, index, orientation);
    let bank = a.radial.as_ref().unwrap();
    a.out = vec![
        a.n as f32,
        a.pulses as f32,
        bank.nf as f32,
        bank.ns as f32,
        a.noise.len() as f32,
        lookup.len() as f32,
        columns.len() as f32,
    ];
    a.out
        .extend(corr.iter().flat_map(|z| [z.r as f32, z.i as f32]));
    a.out
        .extend(columns.iter().map(|&i| f32::from_bits(i as u32)));
    a.out.extend(
        lookup
            .iter()
            .flat_map(|&(i, j)| [f32::from_bits(i as u32), f32::from_bits(j as u32)]),
    );
}

pub(super) fn refine(a: &mut Engine, count: usize) {
    a.radial.as_mut().unwrap().refined = None;
    let mut top: Vec<usize> = (0..a.scores.len())
        .filter(|&i| a.scores[i].is_finite())
        .collect();
    top.sort_unstable_by(|&i, &j| a.scores[j].total_cmp(&a.scores[i]));
    top.truncate(count);
    a.best_score = -1.;
    for i in top {
        let q = values(a, a.radial.as_ref().unwrap().theta(i));
        if q.iter().map(|z| z.power()).sum::<f64>() < 1e-24 {
            continue;
        }
        let mut best = (0., 0.);
        for conj in [false, true] {
            let noise = a.quiet_power(&q, conj);
            if noise > 1e-24 {
                let score = a.event_power(&q, conj) / noise;
                if score > best.0 {
                    best = (score, conj as u8 as f32);
                }
            }
        }

        a.signs[i] = best.1;
        if best.0 > a.best_score {
            a.best_score = best.0;
            a.best = i;
        }
    }
}

fn exact_score(a: &Engine, th: [f64; 3]) -> (f64, f32) {
    if !valid(a, th) {
        return (-1., 0.);
    }
    let q = values(a, th);
    if q.iter().map(|z| z.power()).sum::<f64>() < 1e-24 {
        return (-1., 0.);
    }
    let mut best = (-1., 0.);
    for conj in [false, true] {
        let quiet = a.quiet_power(&q, conj);
        if quiet > 1e-24 {
            let score = a.event_power(&q, conj) / quiet;
            if score.is_finite() && score > best.0 {
                best = (score, conj as u8 as f32);
            }
        }
    }
    best
}
// Bounded Nelder-Mead in grid-spacing units, using the actual quadratic phase.
// The discrete MAX maps remain unchanged; only the reported best-fit is refined.
pub(super) fn refine_peak(a: &mut Engine) {
    if a.best_score < 0. {
        return;
    }
    let bank = a.radial.as_ref().unwrap();
    let seed = bank.theta(a.best);
    let axes = [&bank.axes[0], &bank.axes[2], &bank.axes[1]];
    let bounds = axes.map(|axis| [axis[0], *axis.last().unwrap()]);
    let active: Vec<usize> = (0..3).filter(|&i| axes[i].len() > 1).collect();
    let steps = axes.map(|axis| {
        if axis.len() > 1 {
            axis[1] - axis[0]
        } else {
            1.
        }
    });
    let dimension = active.len();
    if dimension == 0 {
        return;
    }
    let evaluate = |x: Vec<f64>| {
        let mut th = seed;
        for (j, &i) in active.iter().enumerate() {
            th[i] += x[j] * steps[i];
        }
        let value = if (0..3).any(|i| th[i] < bounds[i][0] || th[i] > bounds[i][1]) {
            -1.
        } else {
            exact_score(a, th).0
        };
        (x, value)
    };
    let mut simplex = vec![evaluate(vec![0.; dimension])];
    for (j, &i) in active.iter().enumerate() {
        let mut x = vec![0.; dimension];
        x[j] = if seed[i] + 0.5 * steps[i] <= bounds[i][1] {
            0.5
        } else {
            -0.5
        };
        simplex.push(evaluate(x));
    }
    for _ in 0..300 {
        simplex.sort_by(|a, b| b.1.total_cmp(&a.1));
        let diameter = simplex
            .iter()
            .skip(1)
            .flat_map(|v| v.0.iter().zip(&simplex[0].0).map(|(x, y)| (x - y).abs()))
            .fold(0., f64::max);
        if diameter < 1e-5 {
            break;
        }
        let centre: Vec<f64> = (0..dimension)
            .map(|j| simplex[..dimension].iter().map(|v| v.0[j]).sum::<f64>() / dimension as f64)
            .collect();
        let reflected = evaluate(
            centre
                .iter()
                .zip(&simplex[dimension].0)
                .map(|(c, x)| 2. * c - x)
                .collect(),
        );
        if reflected.1 > simplex[0].1 {
            let expanded = evaluate(
                centre
                    .iter()
                    .zip(&reflected.0)
                    .map(|(c, x)| c + 2. * (x - c))
                    .collect(),
            );
            simplex[dimension] = if expanded.1 > reflected.1 {
                expanded
            } else {
                reflected
            };
        } else if reflected.1 > simplex[dimension - 1].1 {
            simplex[dimension] = reflected;
        } else {
            let outside = reflected.1 > simplex[dimension].1;
            let target = if outside {
                &reflected.0
            } else {
                &simplex[dimension].0
            };
            let contracted = evaluate(
                centre
                    .iter()
                    .zip(target)
                    .map(|(c, x)| c + 0.5 * (x - c))
                    .collect(),
            );
            if contracted.1
                > if outside {
                    reflected.1
                } else {
                    simplex[dimension].1
                }
            {
                simplex[dimension] = contracted;
            } else {
                let best = simplex[0].0.clone();
                for v in simplex.iter_mut().skip(1) {
                    *v = evaluate(
                        v.0.iter()
                            .zip(&best)
                            .map(|(x, b)| b + 0.5 * (x - b))
                            .collect(),
                    );
                }
            }
        }
    }
    simplex.sort_by(|a, b| b.1.total_cmp(&a.1));
    if simplex[0].1 > a.best_score {
        let mut th = seed;
        for (j, &i) in active.iter().enumerate() {
            th[i] += simplex[0].0[j] * steps[i];
        }
        let (score, sign) = exact_score(a, th);
        a.best_score = score;
        a.signs[a.best] = sign;
        a.radial.as_mut().unwrap().refined = Some(th);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn with_engine(f: impl FnOnce(&Engine)) {
        let data = vec![1f32; 1000 * 32 * 2];
        unsafe {
            load(
                data.as_ptr(),
                1000,
                32,
                125,
                12.5e6,
                2e-6,
                25.37e-6,
                77e9,
                1e14,
            );
        }
        assert!(prepare(750, 8, 125, 625, 0) >= 8);
        ENGINE.with(|e| f(e.borrow().as_ref().unwrap()));
    }
    #[test]
    fn corrected_sparse_fft_equals_direct_complex_inner_at_bin() {
        with_engine(|a| {
            let th = [0.3, -320., 5e5];
            let q = values(a, th);
            let nf = 256;
            let ns = 64;
            let f = [a.fs * 17. / nf as f64, 3. / (ns as f64 * a.period)];
            let corr: Vec<Z> = q
                .iter()
                .enumerate()
                .map(|(i, &z)| {
                    let phase = -2.
                        * PI
                        * (f[0] * (i % a.n) as f64 / a.fs + f[1] * (i / a.n) as f64 * a.period);
                    z.mul(Z {
                        r: phase.cos(),
                        i: phase.sin(),
                    })
                })
                .collect();
            for conj in [false, true] {
                let event: Vec<Z> = q
                    .iter()
                    .map(|z| if conj { Z { r: z.r, i: -z.i } } else { *z })
                    .collect();
                let fast = if conj { nf - 17 } else { 17 };
                let slow = if conj { ns - 3 } else { 3 };
                let transformed = spectrum(a, &event, &corr, conj, nf, ns, &[fast]);
                let expected = inner(&q, &event, conj);
                assert!(
                    (transformed[slow].r as f64 - expected.r).abs()
                        < 2e-5 * (1. + expected.power().sqrt())
                );
                assert!(
                    (transformed[slow].i as f64 - expected.i).abs()
                        < 2e-5 * (1. + expected.power().sqrt())
                );
            }
        });
    }
    #[test]
    fn signed_acceleration_full_axes_and_cap() {
        with_engine(|a| {
            let b = [[0.2, 0.4], [-400., 400.], [-1e6, 1e6]];
            let r = Radial::new(a, b, 0.05, 16000000).unwrap();
            assert_eq!([r.axes[0][0], r.axes[0][r.axes[0].len() - 1]], b[0]);
            assert_eq!([r.axes[2][0], r.axes[2][r.axes[2].len() - 1]], b[1]);
            assert_eq!([r.axes[1][0], r.axes[1][r.axes[1].len() - 1]], b[2]);
            assert!(Radial::new(a, b, 0.05, 1).is_none());
            assert!(Radial::new(a, [[0.3, 0.3], [-320., -320.], [0., 0.]], 0.05, 1).is_some());
            assert!(!valid(a, [0.001, 400., 0.]));
        });
    }
}
