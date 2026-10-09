//! Midpoint (range, radial speed, curvature) correction bank for coherent trains.
//! Grid accuracy is sampled and estimated; physical scores and fits are exact.
use super::*;
use std::collections::HashMap;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug)]
struct I(f64, f64);
impl I {
    fn point(x: f64) -> Self {
        Self(x, x)
    }
    fn absmax(self) -> f64 {
        self.0.abs().max(self.1.abs())
    }
    fn sqrt(self) -> Self {
        Self(self.0.max(0.).sqrt(), self.1.sqrt())
    }
    fn reciprocal(self) -> Self {
        assert!(self.0 > 0. || self.1 < 0.);
        Self(1. / self.1, 1. / self.0)
    }
    fn square(self) -> Self {
        Self(
            if self.0 <= 0. && self.1 >= 0. {
                0.
            } else {
                self.0.abs().min(self.1.abs()).powi(2)
            },
            self.absmax().powi(2),
        )
    }
}
impl Add for I {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self(self.0 + b.0, self.1 + b.1)
    }
}
impl Neg for I {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.1, -self.0)
    }
}
impl Sub for I {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + -b
    }
}
impl Mul for I {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        let p = [self.0 * b.0, self.0 * b.1, self.1 * b.0, self.1 * b.1];
        Self(
            p.iter().copied().fold(f64::INFINITY, f64::min),
            p.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    }
}
impl Div for I {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        self * b.reciprocal()
    }
}
#[cfg(test)]
#[derive(Clone, Copy)]
struct D {
    value: I,
    grad: [I; 3],
}
#[cfg(test)]
impl D {
    fn c(x: f64) -> Self {
        Self::interval(I::point(x))
    }
    fn interval(value: I) -> Self {
        Self {
            value,
            grad: [I::point(0.); 3],
        }
    }
    fn variable(value: I, j: usize) -> Self {
        let mut a = Self::interval(value);
        a.grad[j] = I::point(1.);
        a
    }
    fn sqrt(self) -> Self {
        let value = self.value.sqrt();
        Self {
            value,
            grad: self.grad.map(|a| a / (I::point(2.) * value)),
        }
    }
    fn square(self) -> Self {
        Self {
            value: self.value.square(),
            grad: self.grad.map(|a| I::point(2.) * self.value * a),
        }
    }
}
#[cfg(test)]
impl Add for D {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            value: self.value + b.value,
            grad: std::array::from_fn(|j| self.grad[j] + b.grad[j]),
        }
    }
}
#[cfg(test)]
impl Neg for D {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            value: -self.value,
            grad: self.grad.map(|a| -a),
        }
    }
}
#[cfg(test)]
impl Sub for D {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + -b
    }
}
#[cfg(test)]
impl Mul for D {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            value: self.value * b.value,
            grad: std::array::from_fn(|j| self.grad[j] * b.value + self.value * b.grad[j]),
        }
    }
}
#[cfg(test)]
impl Div for D {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        Self {
            value: self.value / b.value,
            grad: std::array::from_fn(|j| {
                (self.grad[j] * b.value - self.value * b.grad[j]) / b.value.square()
            }),
        }
    }
}

#[derive(Clone, Copy)]
struct P {
    value: f64,
    grad: [f64; 3],
}
impl P {
    fn c(value: f64) -> Self {
        Self {
            value,
            grad: [0.; 3],
        }
    }
    fn variable(value: f64, j: usize) -> Self {
        let mut p = Self::c(value);
        p.grad[j] = 1.;
        p
    }
    fn square(self) -> Self {
        Self {
            value: self.value * self.value,
            grad: self.grad.map(|g| 2. * self.value * g),
        }
    }
    fn sqrt(self) -> Self {
        let value = self.value.sqrt();
        Self {
            value,
            grad: self.grad.map(|g| g / (2. * value)),
        }
    }
}
impl Add for P {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            value: self.value + b.value,
            grad: std::array::from_fn(|j| self.grad[j] + b.grad[j]),
        }
    }
}
impl Neg for P {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            value: -self.value,
            grad: self.grad.map(|g| -g),
        }
    }
}
impl Sub for P {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + -b
    }
}
impl Mul for P {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            value: self.value * b.value,
            grad: std::array::from_fn(|j| self.grad[j] * b.value + self.value * b.grad[j]),
        }
    }
}
impl Div for P {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        Self {
            value: self.value / b.value,
            grad: std::array::from_fn(|j| {
                (self.grad[j] * b.value - self.value * b.grad[j]) / (b.value * b.value)
            }),
        }
    }
}
#[derive(Clone)]
pub(super) struct Cell {
    pub bounds: [[f64; 2]; 3], // x, y, v
    pub theta: [f64; 3],
    pub error: f64,
    pub frequencies: [[f64; 2]; 2], // within-ramp Hz, chirp-to-chirp Hz
    pub bank: usize,
    pub candidate: f64,
    pub peak: [f64; 2],
}
struct Bank {
    theta: [f64; 3],
    phase: Vec<f64>,
    amplitude: Vec<f64>,
    frequencies: [[f64; 2]; 2],
}
struct Spectrum {
    event: [HashMap<usize, Z>; 2],
    norm: f64,
}
pub(super) struct Fast {
    banks: Vec<Bank>,
    spectra: HashMap<usize, Spectrum>,
    pending: Vec<([[f64; 2]; 3], f64)>,
    domain: [[f64; 2]; 3],
    pub cells: Vec<Cell>,
    pub limit: f64,
    pub cap: usize,
    pub nf: usize,
    pub ns: usize,
    pub evaluated: usize,
}
#[cfg(test)]
fn range(theta: [D; 3], t: D) -> (D, D, D) {
    let [x, y, v] = theta;
    let beta = v / D::c(C);
    let den = D::c(1.) - beta.square();
    let xx = x - v * t;
    let rr = (xx.square() + den * y.square()).sqrt();
    let r = (rr + beta * xx) / den;
    let r1 = -v * (xx / rr + beta) / den;
    let r2 = v.square() * y.square() / (rr * rr * rr);
    (r, r1, r2)
}
fn point_range(theta: [P; 3], t: P) -> (P, P, P) {
    let [x, y, v] = theta;
    let beta = v / P::c(C);
    let den = P::c(1.) - beta.square();
    let xx = x - v * t;
    let rr = (xx.square() + den * y.square()).sqrt();
    let r = (rr + beta * xx) / den;
    let r1 = -v * (xx / rr + beta) / den;
    let r2 = v.square() * y.square() / (rr * rr * rr);
    (r, r1, r2)
}
fn clocks(a: &Engine) -> (f64, f64, f64) {
    let mid = a.adc + (a.n - 1) as f64 / (2. * a.fs);
    let d = (a.pulses - 1) as f64 * a.period / 2.;
    (mid, mid + d, d + (a.n - 1) as f64 / (2. * a.fs))
}
#[cfg(test)]
fn frequencies(a: &Engine, theta: [D; 3]) -> (D, D, D) {
    let (mid, star, _) = clocks(a);
    let (r, r1, _) = range(theta, D::c(star));
    let slow = -D::c(2. * (a.f0 + a.slope * mid) / C) * r1 + D::c(4. * a.slope / (C * C)) * r * r1;
    let fast = slow - D::c(2. * a.slope / C) * r;
    let phi0 = -D::c(4. * PI * (a.f0 + a.slope * mid) / C) * r
        + D::c(4. * PI * a.slope / (C * C)) * r.square();
    (fast, slow, phi0)
}
fn point_frequencies(a: &Engine, theta: [P; 3]) -> (P, P, P) {
    let (mid, star, _) = clocks(a);
    let (r, r1, _) = point_range(theta, P::c(star));
    let slow = -P::c(2. * (a.f0 + a.slope * mid) / C) * r1 + P::c(4. * a.slope / (C * C)) * r * r1;
    let fast = slow - P::c(2. * a.slope / C) * r;
    let phi0 = -P::c(4. * PI * (a.f0 + a.slope * mid) / C) * r
        + P::c(4. * PI * a.slope / (C * C)) * r.square();
    (fast, slow, phi0)
}
fn scalar_frequencies(a: &Engine, theta: [f64; 3]) -> [f64; 3] {
    let (f, b, p) = point_frequencies(a, theta.map(P::c));
    [f.value, b.value, p.value]
}
#[cfg(test)]
fn conservative_cell_bound(a: &Engine, bounds: [[f64; 2]; 3]) -> (f64, [f64; 3], [[f64; 2]; 2]) {
    let th = std::array::from_fn(|j| D::variable(I(bounds[j][0], bounds[j][1]), j));
    let (mid, star, h) = clocks(a);
    let reset = (a.pulses - 1) as f64 * a.period / 2.;
    let end = (a.pulses - 1) as f64 * a.period + a.adc + (a.n - 1) as f64 / a.fs;
    let mut b2 = [0f64; 3];
    let mut b1 = [0f64; 3];
    for k in 0..16 {
        let t = D::interval(I(
            a.adc + (end - a.adc) * k as f64 / 16.,
            a.adc + (end - a.adc) * (k + 1) as f64 / 16.,
        ));
        let (r, r1, r2) = range(th, t);
        let w = D::c(2. * PI) * (D::c(a.f0) + D::c(a.slope) * (t - D::c(star) + D::c(mid)));
        let second = -D::c(2. / C) * (w * r2 + D::c(4. * PI * a.slope) * r1)
            + D::c(8. * PI * a.slope / (C * C)) * (r1.square() + r * r2);
        for j in 0..3 {
            b2[j] = b2[j].max(second.grad[j].absmax());
            b1[j] = b1[j].max(r1.grad[j].absmax());
        }
    }
    let parts = std::array::from_fn(|j| {
        (bounds[j][1] - bounds[j][0]) / 2.
            * (h * h / 2. * b2[j] + 4. * PI * a.slope.abs() / C * reset * h * b1[j])
            * (1. + 1e-10)
    });
    let (f, s, _) = frequencies(a, th);
    (
        parts.iter().sum(),
        parts,
        [[f.value.0, f.value.1], [s.value.0, s.value.1]],
    )
}
// Range, radial velocity and radial acceleration at the train midpoint.
fn radial_theta(a: &Engine, c: [P; 3], sign: f64) -> [P; 3] {
    let [r, u, g] = c;
    let v = (u.square() + g * r).sqrt() * P::c(sign);
    let (_, star, _) = clocks(a);
    [
        v * P::c(star) - u * r / v,
        (g * r * r * r / v.square()).sqrt(),
        v,
    ]
}
fn radial_coords(a: &Engine, th: [f64; 3]) -> [f64; 3] {
    let (_, star, _) = clocks(a);
    let [x, y, v] = th;
    let xx = v * star - x;
    let r = xx.hypot(y);
    [r, v * xx / r, v * v * y * y / (r * r * r)]
}
fn inside(th: [f64; 3], domain: [[f64; 2]; 3]) -> bool {
    (0..3).all(|j| {
        th[j].is_finite() && th[j] >= domain[j][0] - 1e-12 && th[j] <= domain[j][1] + 1e-12
    })
}
fn physical_bounds(
    a: &Engine,
    b: [[f64; 2]; 3],
    sign: f64,
    domain: [[f64; 2]; 3],
) -> Option<[[f64; 2]; 3]> {
    let [r, u, g] = b.map(|b| I(b[0], b[1]));
    let den = u.square() + g * r;
    let mut v = den.sqrt();
    let vbound = if sign > 0. {
        I(domain[2][0].max(0.), domain[2][1])
    } else {
        I((-domain[2][1]).max(0.), -domain[2][0])
    };
    v = I(v.0.max(vbound.0), v.1.min(vbound.1));
    if v.0 > v.1 || v.1 <= 0. {
        return None;
    }
    v.0 = v.0.max(1e-10);
    let (_, star, _) = clocks(a);
    let x = I::point(sign) * (v * I::point(star) - u * r / v);
    let y = (g * r * r * r / den).sqrt();
    let v = I::point(sign) * v;
    let result = [x, y, v].map(|q| [q.0, q.1]);
    let result = std::array::from_fn(|j| {
        [
            result[j][0].max(domain[j][0]),
            result[j][1].min(domain[j][1]),
        ]
    });
    if result.iter().any(|b| b[0] > b[1]) {
        None
    } else {
        Some(result)
    }
}
// A sampled local metric, with a safety factor; not a global interval certificate.
fn cell_bound(
    a: &Engine,
    bounds: [[f64; 2]; 3],
    sign: f64,
    domain: [[f64; 2]; 3],
    hull: [[f64; 2]; 3],
) -> Option<(f64, [f64; 3], [[f64; 2]; 2], [f64; 3])> {
    let centre = bounds.map(|b| (b[0] + b[1]) / 2.);
    let mut points: Vec<[f64; 3]> = vec![];
    for x in 0..3 {
        for y in 0..3 {
            for v in 0..3 {
                let c = std::array::from_fn(|j| match [x, y, v][j] {
                    0 => bounds[j][0],
                    1 => centre[j],
                    _ => bounds[j][1],
                });
                if inside(radial_theta(a, c.map(P::c), sign).map(|p| p.value), domain) {
                    points.push(c);
                }
            }
        }
    }
    // Physical hull and original-domain nodes help with cells on curved boundaries.
    for box_ in [hull, domain] {
        for bits in 0..9 {
            let th = if bits == 8 {
                box_.map(|b| (b[0] + b[1]) / 2.)
            } else {
                std::array::from_fn(|j| box_[j][(bits >> j) & 1])
            };
            if th[2] * sign <= 0. {
                continue;
            }
            let c = radial_coords(a, th);
            if (0..3).all(|j| c[j] >= bounds[j][0] && c[j] <= bounds[j][1]) {
                points.push(c);
            }
        }
    }
    if points.is_empty() {
        for start in [centre, radial_coords(a, hull.map(|b| (b[0] + b[1]) / 2.))] {
            let mut c = start;
            for _ in 0..64 {
                for j in 0..3 {
                    c[j] = c[j].clamp(bounds[j][0], bounds[j][1]);
                }
                let th = radial_theta(a, c.map(P::c), sign).map(|p| p.value);
                if inside(th, domain) {
                    points.push(c);
                    break;
                }
                let projected = std::array::from_fn(|j| th[j].clamp(domain[j][0], domain[j][1]));
                c = radial_coords(a, projected);
            }
        }
        if points.is_empty() {
            return None;
        }
    }
    let representative = points
        .iter()
        .min_by(|a, b| {
            let d = |c: &[f64; 3]| {
                (0..3)
                    .map(|j| {
                        ((c[j] - centre[j]) / (bounds[j][1] - bounds[j][0]).max(1e-30)).powi(2)
                    })
                    .sum::<f64>()
            };
            d(a).total_cmp(&d(b))
        })
        .unwrap();
    let theta = radial_theta(a, representative.map(P::c), sign).map(|p| p.value);
    let centre_curves = curves(a, theta).0;
    let (mid, _, _) = clocks(a);
    let mut error = 0f64;
    let mut derivatives = [0f64; 3];
    let mut freq = [[f64::INFINITY, f64::NEG_INFINITY]; 2];
    for point in &points {
        let c = std::array::from_fn(|j| P::variable(point[j], j));
        let th = radial_theta(a, c, sign);
        let (fast, slow, p0) = point_frequencies(a, th);
        for (k, f) in [fast, slow].iter().enumerate() {
            freq[k][0] = freq[k][0].min(f.value);
            freq[k][1] = freq[k][1].max(f.value);
        }
        for k in 0..a.pulses {
            for j in [0, (a.n - 1) / 4, (a.n - 1) / 2, 3 * (a.n - 1) / 4, a.n - 1] {
                let u = a.adc + j as f64 / a.fs;
                let d = (k as f64 - (a.pulses - 1) as f64 / 2.) * a.period;
                let (r, _, _) = point_range(th, P::c(k as f64 * a.period + u));
                let phase = -P::c(4. * PI * (a.f0 + a.slope * u) / C) * r
                    + P::c(4. * PI * a.slope / (C * C)) * r.square();
                let residual = phase - p0 - P::c(2. * PI) * (fast * P::c(u - mid) + slow * P::c(d));
                error = error.max((residual.value - centre_curves[k * a.n + j]).abs());
                for p in 0..3 {
                    derivatives[p] = derivatives[p].max(residual.grad[p].abs());
                }
            }
        }
    }
    let parts = std::array::from_fn(|j| 1.25 * (bounds[j][1] - bounds[j][0]) / 2. * derivatives[j]);
    Some((1.25 * error, parts, freq, theta))
}
impl Fast {
    pub fn new(a: &Engine, bounds: [[f64; 2]; 3], loss: f64, cap: usize) -> Option<Self> {
        if !loss.is_finite()
            || !(0. < loss && loss < 1.)
            || cap == 0
            || bounds
                .iter()
                .any(|b| !b[0].is_finite() || !b[1].is_finite() || b[0] > b[1])
            || bounds[1][0] <= 0.
            || bounds[2].iter().any(|v| v.abs() > 100000.)
        {
            return None;
        }
        Some(Self {
            banks: vec![],
            spectra: HashMap::new(),
            pending: {
                let (_, star, _) = clocks(a);
                let vmax = bounds[2][0].abs().max(bounds[2][1].abs());
                let xmax = (bounds[0][0] - bounds[2][1] * star)
                    .abs()
                    .max((bounds[0][1] - bounds[2][0] * star).abs());
                let rmax = xmax.hypot(bounds[1][1]);
                let vmin = if bounds[2][0] <= 0. && bounds[2][1] >= 0. {
                    0.
                } else {
                    bounds[2][0].abs().min(bounds[2][1].abs())
                };
                let radial = [
                    [bounds[1][0], rmax],
                    [-vmax, vmax],
                    [
                        (vmin * vmin * bounds[1][0] * bounds[1][0] / rmax.powi(3)).max(1e-8),
                        vmax * vmax / bounds[1][0] + 1e-8,
                    ],
                ];
                let mut p = vec![];
                if bounds[2][1] > 0. {
                    p.push((radial, 1.));
                }
                if bounds[2][0] < 0. {
                    p.push((radial, -1.));
                }
                p
            },
            domain: bounds,
            cells: if bounds[2][0] <= 0. && bounds[2][1] >= 0. {
                let theta = [
                    (bounds[0][0] + bounds[0][1]) / 2.,
                    (bounds[1][0] + bounds[1][1]) / 2.,
                    0.,
                ];
                let xmin = if bounds[0][0] <= 0. && bounds[0][1] >= 0. {
                    0.
                } else {
                    bounds[0][0].abs().min(bounds[0][1].abs())
                };
                let xmax = bounds[0][0].abs().max(bounds[0][1].abs());
                let rmin = xmin.hypot(bounds[1][0]);
                let rmax = xmax.hypot(bounds[1][1]);
                let mut b = bounds;
                b[2] = [0., 0.];
                vec![Cell {
                    bounds: b,
                    theta,
                    error: 0.,
                    frequencies: [
                        [-2. * a.slope * rmax / C, -2. * a.slope * rmin / C],
                        [0., 0.],
                    ],
                    bank: 0,
                    candidate: 0.,
                    peak: [0.; 2],
                }]
            } else {
                vec![]
            },
            limit: (1. - loss).sqrt().acos(),
            cap,
            nf: (a.n * 4).next_power_of_two(),
            ns: if a.pulses == 1 {
                1
            } else {
                (a.pulses * 4).next_power_of_two()
            },
            evaluated: 0,
        })
    }
    pub fn build(&mut self, a: &Engine, count: usize) -> isize {
        for _ in 0..count {
            let Some((bounds, sign)) = self.pending.pop() else {
                self.group(a);
                return if self.cells.is_empty() {
                    -3
                } else {
                    self.cells.len() as isize
                };
            };
            let Some(hull) = physical_bounds(a, bounds, sign, self.domain) else {
                continue;
            };
            let result = cell_bound(a, bounds, sign, self.domain, hull);
            let Some((error, parts, freq, theta)) = result else {
                continue;
            };

            if error <= self.limit {
                self.cells.push(Cell {
                    bounds: hull,
                    theta,
                    error,
                    frequencies: freq,
                    bank: 0,
                    candidate: 0.,
                    peak: [0.; 2],
                });
            } else {
                let j = (0..3)
                    .max_by(|&i, &j| parts[i].total_cmp(&parts[j]))
                    .unwrap();
                let mid = (bounds[j][0] + bounds[j][1]) / 2.;
                if mid == bounds[j][0] || mid == bounds[j][1] {
                    continue;
                }
                let (mut left, mut right) = (bounds, bounds);
                left[j][1] = mid;
                right[j][0] = mid;
                self.pending.push((right, sign));
                self.pending.push((left, sign));
            }
            if self.cells.len() + self.pending.len() > self.cap {
                return -1;
            }
        }
        if self.pending.is_empty() {
            self.group(a);
            if self.cells.is_empty() {
                -3
            } else {
                self.cells.len() as isize
            }
        } else {
            0
        }
    }
    pub fn bank_count(&self) -> usize {
        self.banks.len()
    }
    fn group(&mut self, a: &Engine) {
        if !self.banks.is_empty() {
            return;
        }
        let mut keys: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        for cell in &mut self.cells {
            let (phase, amplitude) = curves(a, cell.theta);
            let signature = [phase[0], phase[a.n - 1], *phase.last().unwrap()];
            let key = signature.map(|v| (v / (2. * self.limit)).floor() as i64);
            let allowance = self.limit;
            let mut found = None;
            'neighbours: for x in -1..=1 {
                for y in -1..=1 {
                    for z in -1..=1 {
                        let k = [key[0] + x, key[1] + y, key[2] + z];
                        if let Some(indices) = keys.get(&k) {
                            for &index in indices {
                                let bank = &self.banks[index];
                                if phase
                                    .iter()
                                    .zip(&bank.phase)
                                    .any(|(a, b)| (a - b).abs() > allowance)
                                {
                                    continue;
                                }
                                if amplitude
                                    .iter()
                                    .zip(&bank.amplitude)
                                    .any(|(a, b)| (a / b - 1.).abs() > 0.025)
                                {
                                    continue;
                                }
                                found = Some(index);
                                break 'neighbours;
                            }
                        }
                    }
                }
            }
            cell.bank = found.unwrap_or_else(|| {
                let index = self.banks.len();
                self.banks.push(Bank {
                    theta: cell.theta,
                    phase,
                    amplitude,
                    frequencies: cell.frequencies,
                });
                keys.entry(key).or_default().push(index);
                index
            });
            for j in 0..2 {
                self.banks[cell.bank].frequencies[j][0] =
                    self.banks[cell.bank].frequencies[j][0].min(cell.frequencies[j][0]);
                self.banks[cell.bank].frequencies[j][1] =
                    self.banks[cell.bank].frequencies[j][1].max(cell.frequencies[j][1]);
            }
        }
    }
    pub fn pending(&self) -> usize {
        self.pending.len()
    }
}
fn curves(a: &Engine, theta: [f64; 3]) -> (Vec<f64>, Vec<f64>) {
    let (mid, star, _) = clocks(a);
    let f = scalar_frequencies(a, theta);
    let th = theta.map(P::c);
    let rstar = point_range(th, P::c(star)).0.value;
    let mut phase = Vec::with_capacity(a.n * a.pulses);
    let mut amp = Vec::with_capacity(a.n * a.pulses);
    for k in 0..a.pulses {
        for j in 0..a.n {
            let u = a.adc + j as f64 / a.fs;
            let d = (k as f64 - (a.pulses - 1) as f64 / 2.) * a.period;
            let r = point_range(th, P::c(k as f64 * a.period + u)).0.value;
            phase.push(
                -4. * PI * (a.f0 + a.slope * u) * r / C + 4. * PI * a.slope * r * r / (C * C)
                    - f[2]
                    - 2. * PI * (f[0] * (u - mid) + f[1] * d),
            );
            amp.push((rstar / r).powi(2));
        }
    }
    (phase, amp)
}
fn correction(a: &Engine, theta: [f64; 3]) -> Vec<Z> {
    let q: Vec<Z> = template(
        a.n, a.pulses, a.fs, a.adc, a.period, a.f0, a.slope, theta[0], theta[1], theta[2],
        a.receiver,
    )
    .into_iter()
    .map(|(r, i)| Z { r, i })
    .collect();
    let f = scalar_frequencies(a, theta);
    let (mid, _, _) = clocks(a);
    q.iter()
        .enumerate()
        .map(|(i, &z)| {
            let eta = a.adc + (i % a.n) as f64 / a.fs - mid;
            let d = ((i / a.n) as f64 - (a.pulses - 1) as f64 / 2.) * a.period;
            let p = f[2] + 2. * PI * (f[0] * eta + f[1] * d);
            z.mul(Z {
                r: p.cos(),
                i: -p.sin(),
            })
        })
        .collect()
}
fn transform(a: &Engine, z: &[Z], correction: &[Z], conj: bool, nf: usize, ns: usize) -> Vec<Z> {
    let mut out = vec![Z::default(); nf * ns];
    for k in 0..a.pulses {
        for j in 0..a.n {
            let q = correction[k * a.n + j];
            let q = Z {
                r: q.r,
                i: if conj { q.i } else { -q.i },
            };
            out[k * nf + j] = z[k * a.n + j].mul(q);
        }
        fft(&mut out[k * nf..(k + 1) * nf]);
    }
    if ns > 1 {
        let mut column = vec![Z::default(); ns];
        for j in 0..nf {
            for k in 0..ns {
                column[k] = out[k * nf + j];
            }
            fft(&mut column);
            for k in 0..ns {
                out[k * nf + j] = column[k];
            }
        }
    }
    out
}
fn lift_frequency(f: f64, period: f64, interval: [f64; 2], half_bin: f64) -> Option<f64> {
    let f = f + period * ((interval[0] + interval[1]) / 2. / period - f / period).round();
    if f >= interval[0] - half_bin && f <= interval[1] + half_bin {
        Some(f)
    } else {
        None
    }
}
pub(super) fn physical_score(a: &Engine, theta: [f64; 3]) -> (f64, f32) {
    let q: Vec<Z> = template(
        a.n, a.pulses, a.fs, a.adc, a.period, a.f0, a.slope, theta[0], theta[1], theta[2],
        a.receiver,
    )
    .into_iter()
    .map(|(r, i)| Z { r, i })
    .collect();
    score_q(a, &q)
}
fn score_q(a: &Engine, q: &[Z]) -> (f64, f32) {
    let mut best = (0., 0.);
    for conj in [false, true] {
        let denom = a
            .noise
            .iter()
            .map(|z| inner(q, z, conj).power())
            .sum::<f64>()
            / a.noise.len() as f64;
        if denom > 1e-24 {
            let value = inner(q, &a.event, conj).power() / denom;
            if value > best.0 {
                best = (value, conj as u8 as f32);
            }
        }
    }
    best
}
pub(super) fn search_cell(a: &mut Engine, index: usize) {
    let fast = a.fast.as_ref().unwrap();
    let cell = fast.cells[index].clone();
    let (nf, ns) = (fast.nf, fast.ns);
    let q = a.q(index);
    let freq = scalar_frequencies(a, cell.theta);
    if !a.fast.as_ref().unwrap().spectra.contains_key(&cell.bank) {
        let theta = a.fast.as_ref().unwrap().banks[cell.bank].theta;
        let c = correction(a, theta);
        let intervals = a.fast.as_ref().unwrap().banks[cell.bank].frequencies;
        let mut spectrum = Spectrum {
            event: std::array::from_fn(|_| HashMap::new()),
            norm: c.iter().map(|z| z.power()).sum(),
        };
        for orientation in 0..2 {
            let conj = orientation == 1;
            let intervals = intervals.map(|b| if conj { [-b[1], -b[0]] } else { b });
            let event = transform(a, &a.event, &c, conj, nf, ns);
            for k in 0..ns {
                for j in 0..nf {
                    let f = (if j < nf / 2 {
                        j as f64
                    } else {
                        j as f64 - nf as f64
                    }) * a.fs
                        / nf as f64;
                    let slow = (if k < ns / 2 {
                        k as f64
                    } else {
                        k as f64 - ns as f64
                    }) / (ns as f64 * a.period);
                    if lift_frequency(f, a.fs, intervals[0], a.fs / (2. * nf as f64)).is_some()
                        && (ns == 1
                            || lift_frequency(
                                slow,
                                1. / a.period,
                                intervals[1],
                                1. / (2. * ns as f64 * a.period),
                            )
                            .is_some())
                    {
                        spectrum.event[orientation].insert(k * nf + j, event[k * nf + j]);
                    }
                }
            }
        }
        a.fast.as_mut().unwrap().spectra.insert(cell.bank, spectrum);
    }
    let (mut candidate, mut peak) = (0., [freq[0], freq[1]]);
    for conj in [false, true] {
        let sign = if conj { -1. } else { 1. };
        let intervals = cell
            .frequencies
            .map(|b| if conj { [-b[1], -b[0]] } else { b });
        let spectrum = &a.fast.as_ref().unwrap().spectra[&cell.bank];
        let orientation = conj as usize;
        let event = &spectrum.event[orientation];
        for k in 0..ns {
            let slow = (if k < ns / 2 {
                k as f64
            } else {
                k as f64 - ns as f64
            }) / (ns as f64 * a.period);
            let slow = if ns == 1 {
                Some(sign * freq[1])
            } else {
                lift_frequency(
                    slow,
                    1. / a.period,
                    intervals[1],
                    1. / (2. * ns as f64 * a.period),
                )
            };
            let Some(slow) = slow else {
                continue;
            };
            for j in 0..nf {
                let f = (if j < nf / 2 {
                    j as f64
                } else {
                    j as f64 - nf as f64
                }) * a.fs
                    / nf as f64;
                let Some(f) = lift_frequency(f, a.fs, intervals[0], a.fs / (2. * nf as f64)) else {
                    continue;
                };
                let id = k * nf + j;
                if let Some(z) = event.get(&id) {
                    let value = z.power() / spectrum.norm;
                    if value > candidate {
                        candidate = value;
                        peak = [sign * f, sign * slow];
                    }
                }
            }
        }
    }
    let (value, sign) = score_q(a, &q);
    a.scores[index] = value as f32;
    a.signs[index] = sign;
    if value > a.best_score {
        a.best_score = value;
        a.best = index;
    }
    // Time scans retain the original exact physical-template statistic.
    for conj in [false, true] {
        let denom = a
            .noise
            .iter()
            .map(|z| inner(&q, z, conj).power())
            .sum::<f64>()
            / a.noise.len() as f64;
        if denom > 1e-24 {
            for (j, z) in a.scan_data.iter().enumerate() {
                let value = (inner(&q, z, conj).power() / denom) as f32;
                if value > a.scan_peaks[j] {
                    a.scan_peaks[j] = value;
                    a.scan_theta[j] = index;
                }
            }
        }
    }
    let fast = a.fast.as_mut().unwrap();
    fast.cells[index].candidate = candidate;
    fast.cells[index].peak = peak;
    fast.evaluated += 1;
}
fn physical_aliases(a: &Engine, cell: &Cell, peak: [f64; 2]) -> Vec<[f64; 2]> {
    let sf = 1. / a.period;
    let half = sf / (2. * a.fast.as_ref().unwrap().ns as f64);
    let lo = ((cell.frequencies[1][0] - half - peak[1]) / sf).ceil() as i64;
    let hi = ((cell.frequencies[1][1] + half - peak[1]) / sf).floor() as i64;
    let (_, star, _) = clocks(a);
    let xb = cell.bounds[0];
    let vb = cell.bounds[2];
    let yb = cell.bounds[1];
    let rmax = (xb[0] - vb[1] * star)
        .abs()
        .max((xb[1] - vb[0] * star).abs())
        .hypot(yb[1]);
    let vmax = vb[0].abs().max(vb[1].abs());
    let mut out = vec![];
    for k in lo..=hi {
        let slow = peak[1] + k as f64 * sf;
        if a.slope.abs() > 1e-20 {
            let r = C * (slow - peak[0]) / (2. * a.slope);
            let radial =
                slow / (-2. * (a.f0 + a.slope * clocks(a).0) / C + 4. * a.slope * r / (C * C));
            let margin = 1e-4 + rmax * vmax / C * 2.;
            if r < yb[0] - margin
                || r > rmax + margin
                || radial.abs() > vmax * (1. + 2. * vmax / C) + 1.
            {
                continue;
            }
        }
        out.push([peak[0], slow]);
    }
    if out.is_empty() {
        out.push(peak);
    }
    out
}
pub(super) fn refine_cell(a: &mut Engine, index: usize) {
    let cell = a.fast.as_ref().unwrap().cells[index].clone();
    let aliases = physical_aliases(a, &cell, cell.peak);
    for peak in aliases {
        a.fast.as_mut().unwrap().cells[index].peak = peak;
        refine_one(a, index);
    }
}
fn refine_one(a: &mut Engine, index: usize) {
    let cell = a.fast.as_ref().unwrap().cells[index].clone();
    let widths = cell.bounds.map(|b| b[1] - b[0]);
    let mut theta = cell.theta;
    // Bring the physical frequency pair towards the FFT peak, inside this cell.
    for _ in 0..8 {
        let f = scalar_frequencies(a, theta);
        let mut jac = [[0f64; 3]; 2];
        for j in 0..3 {
            if widths[j] == 0. {
                continue;
            }
            let mut th = theta;
            th[j] += widths[j] * 1e-4;
            let f1 = scalar_frequencies(a, th);
            for k in 0..2 {
                jac[k][j] = (f1[k] - f[k]) / 1e-4;
            }
        }
        let residual = [
            cell.peak[0] - f[0],
            if a.pulses == 1 {
                0.
            } else {
                cell.peak[1] - f[1]
            },
        ];
        let rows = if a.pulses == 1 || jac[1].iter().map(|x| x * x).sum::<f64>() < 1e-20 {
            1
        } else {
            2
        };
        let aa = jac[0].iter().map(|x| x * x).sum::<f64>();
        let bb = jac[1].iter().map(|x| x * x).sum::<f64>();
        let ab = (0..3).map(|j| jac[0][j] * jac[1][j]).sum::<f64>();
        let (w0, w1) = if rows == 1 {
            if aa < 1e-20 {
                break;
            }
            (residual[0] / aa, 0.)
        } else {
            let determinant = aa * bb - ab * ab;
            if determinant <= 1e-12 * aa * bb {
                break;
            }
            (
                (bb * residual[0] - ab * residual[1]) / determinant,
                (aa * residual[1] - ab * residual[0]) / determinant,
            )
        };
        for j in 0..3 {
            theta[j] = (theta[j] + widths[j] * (jac[0][j] * w0 + jac[1][j] * w1).clamp(-0.5, 0.5))
                .clamp(cell.bounds[j][0], cell.bounds[j][1]);
        }
    }
    let (mut value, mut orientation) = physical_score(a, theta);
    if a.scores[index] as f64 > value {
        theta = cell.theta;
        value = a.scores[index] as f64;
        orientation = a.signs[index];
    }
    for round in 0..10 {
        let step = 0.25 * 0.5f64.powi(round / 2);
        for j in 0..3 {
            for direction in [-1., 1.] {
                let mut trial = theta;
                trial[j] = (trial[j] + direction * step * widths[j])
                    .clamp(cell.bounds[j][0], cell.bounds[j][1]);
                if trial == theta {
                    continue;
                }
                let (score, sign) = physical_score(a, trial);
                if score > value {
                    theta = trial;
                    value = score;
                    orientation = sign;
                }
            }
        }
    }
    a.fast.as_mut().unwrap().cells[index].theta = theta;
    a.scores[index] = value as f32;
    a.signs[index] = orientation;
    if value > a.best_score {
        a.best = index;
        a.best_score = value;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn with_clock(pulses: usize, f: impl FnOnce(&Engine)) {
        let data = vec![1f32; 1000 * 128 * 2];
        unsafe {
            load(
                data.as_ptr(),
                1000,
                128,
                125,
                12.5e6,
                2e-6,
                15.71e-6,
                79.5e9,
                60e12,
            );
        }
        assert!(prepare(750, pulses, 125, 625, 0) >= 8);
        ENGINE.with(|e| f(e.borrow().as_ref().unwrap()));
    }
    fn phase_correction(a: &Engine, theta: [f64; 3], k: usize, j: usize) -> f64 {
        let u = a.adc + j as f64 / a.fs;
        let t = k as f64 * a.period + u;
        let (r, _, _) = range(theta.map(D::c), D::c(t));
        let r = r.value.0;
        let phase = -4. * PI * (a.f0 + a.slope * u) * r / C + 4. * PI * a.slope * r * r / (C * C);
        let (mid, _, _) = clocks(a);
        let f = scalar_frequencies(a, theta);
        phase
            - f[2]
            - 2. * PI
                * (f[0] * (u - mid) + f[1] * (k as f64 - (a.pulses - 1) as f64 / 2.) * a.period)
    }
    #[test]
    fn midpoint_inverse_and_derivatives() {
        with_clock(16, |a| {
            for th in [
                [0.3, 0.1, 341.51],
                [-0.3, 0.1, -341.51],
                [0.05, 0.2, -120.],
                [0.001, 0.05, 500.],
            ] {
                let c = radial_coords(a, th);
                let sign = th[2].signum();
                let back = radial_theta(a, c.map(P::c), sign);
                for j in 0..3 {
                    assert!((back[j].value - th[j]).abs() < 1e-10);
                }
                let (_, star, _) = clocks(a);
                for h in [-50e-6, 0., 50e-6] {
                    let actual = (th[2] * (star + h) - th[0]).hypot(th[1]);
                    let r =
                        (c[0] * c[0] + 2. * c[0] * c[1] * h + (c[1] * c[1] + c[0] * c[2]) * h * h)
                            .sqrt();
                    assert!((r - actual).abs() < 1e-12);
                }
                let derivative =
                    radial_theta(a, std::array::from_fn(|j| P::variable(c[j], j)), sign);
                for j in 0..3 {
                    let step = c[j].abs().max(1.) * 1e-5;
                    let mut left = c;
                    let mut right = c;
                    left[j] -= step;
                    right[j] += step;
                    let left = radial_theta(a, left.map(P::c), sign);
                    let right = radial_theta(a, right.map(P::c), sign);
                    for p in 0..3 {
                        let numerical = (right[p].value - left[p].value) / (2. * step);
                        assert!(
                            (numerical - derivative[p].grad[j]).abs()
                                < 1e-6 * numerical.abs().max(1.)
                        );
                    }
                }
            }
        });
    }
    #[test]
    fn continuous_bounds_include_retardation_and_resets() {
        for pulses in [1, 2, 4, 8, 16] {
            with_clock(pulses, |a| {
                for bounds in [
                    [[0.29, 0.31], [0.09, 0.11], [330., 350.]],
                    [[0.29, 0.31], [0.09, 0.11], [-7000., -6900.]],
                ] {
                    let (bound, _, _) = conservative_cell_bound(a, bounds);
                    let centre = bounds.map(|b| (b[0] + b[1]) / 2.);
                    for bits in 0..8 {
                        let th = std::array::from_fn(|j| bounds[j][(bits >> j) & 1]);
                        for k in 0..pulses {
                            for j in [0, 32, 64, 127] {
                                let error = (phase_correction(a, th, k, j)
                                    - phase_correction(a, centre, k, j))
                                .abs();
                                assert!(error <= bound + 1e-8, "{pulses}: {error} > {bound}");
                            }
                        }
                    }
                }
            });
        }
    }
    #[test]
    fn two_dimensional_fft_matches_direct_physical_template() {
        with_clock(4, |a| {
            let th = [0.3, 0.1, 341.51];
            let q: Vec<Z> = template(
                a.n, a.pulses, a.fs, a.adc, a.period, a.f0, a.slope, th[0], th[1], th[2], false,
            )
            .into_iter()
            .map(|(r, i)| Z { r, i })
            .collect();
            let (nf, ns) = (512, 16);
            let f = scalar_frequencies(a, th);
            let (mid, _, _) = clocks(a);
            let correction: Vec<Z> = q
                .iter()
                .enumerate()
                .map(|(i, &z)| {
                    let eta = a.adc + (i % a.n) as f64 / a.fs - mid;
                    let d = ((i / a.n) as f64 - 1.5) * a.period;
                    let phase = f[2] + 2. * PI * (f[0] * eta + f[1] * d);
                    z.mul(Z {
                        r: phase.cos(),
                        i: -phase.sin(),
                    })
                })
                .collect();
            for conj in [false, true] {
                let out = transform(a, &a.event, &correction, conj, nf, ns);
                for (k, j) in [(0, 0), (3, 37), (15, 500)] {
                    let f = j as f64 * a.fs / nf as f64;
                    let slow = k as f64 / (ns as f64 * a.period);
                    let qbin: Vec<Z> = correction
                        .iter()
                        .enumerate()
                        .map(|(i, &c)| {
                            let phase = 2.
                                * PI
                                * (f * (i % a.n) as f64 / a.fs
                                    + slow * (i / a.n) as f64 * a.period);
                            let c = if conj { Z { r: c.r, i: -c.i } } else { c };
                            c.mul(Z {
                                r: phase.cos(),
                                i: phase.sin(),
                            })
                        })
                        .collect();
                    let direct = inner(&qbin, &a.event, false);
                    assert!((out[k * nf + j].r - direct.r).abs() < 1e-7);
                    assert!((out[k * nf + j].i - direct.i).abs() < 1e-7);
                }
            }
        });
    }
}
