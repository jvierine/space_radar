//! Browser DSP: complex quiet mean, chirp FFT and retarded Eq.25 filter bank.
//! Memory crossing the JS ABI is float32; phase, sums and ratios use float64.
use std::cell::RefCell;
use std::f64::consts::PI;
const C: f64 = 299_792_458.0;
#[derive(Clone, Copy, Default)]
struct Z {
    r: f64,
    i: f64,
}
impl Z {
    fn power(self) -> f64 {
        self.r * self.r + self.i * self.i
    }
    fn mul(self, b: Self) -> Self {
        Self {
            r: self.r * b.r - self.i * b.i,
            i: self.r * b.i + self.i * b.r,
        }
    }
    fn hp(f: f64, corner: f64) -> Self {
        let d = corner * corner + f * f;
        Self {
            r: f * f / d,
            i: f * corner / d,
        }
    }
}
fn fft(a: &mut [Z]) {
    let n = a.len();
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let w = Z {
            r: (-2.0 * PI / len as f64).cos(),
            i: (-2.0 * PI / len as f64).sin(),
        };
        for k in (0..n).step_by(len) {
            let mut wi = Z { r: 1.0, i: 0.0 };
            for t in 0..len / 2 {
                let u = a[k + t];
                let v = a[k + t + len / 2].mul(wi);
                a[k + t] = Z {
                    r: u.r + v.r,
                    i: u.i + v.i,
                };
                a[k + t + len / 2] = Z {
                    r: u.r - v.r,
                    i: u.i - v.i,
                };
                wi = wi.mul(w);
            }
        }
        len *= 2;
    }
}
fn inner(q: &[Z], z: &[Z], conjugated: bool) -> Z {
    let mut r = 0.;
    let mut i = 0.;
    for (a, b) in q.iter().zip(z) {
        let qi = if conjugated { -a.i } else { a.i };
        r += a.r * b.r + qi * b.i;
        i += a.r * b.i - qi * b.r;
    }
    Z { r, i }
}
pub fn template(
    samples: usize,
    pulses: usize,
    fs: f64,
    adc: f64,
    period: f64,
    f0: f64,
    slope: f64,
    x0: f64,
    y0: f64,
    v: f64,
    receiver: bool,
) -> Vec<(f64, f64)> {
    let mut q = Vec::with_capacity(samples * pulses);
    for k in 0..pulses {
        for j in 0..samples {
            let u = adc + j as f64 / fs;
            let t = k as f64 * period + u;
            let x = v * t - x0;
            let rr = x * x + y0 * y0;
            let d = rr / ((x * v * x * v + (C * C - v * v) * rr).sqrt() + x * v);
            let tau = 2. * d;
            let xs = x - v * d;
            let r = xs.hypot(y0);
            let vr = v * xs / r;
            let td = 2. * vr / (C + vr);
            let phase = 2. * PI * (-f0 * tau - slope * u * tau + 0.5 * slope * tau * tau);
            let f = -slope * tau - (f0 + slope * (u - tau)) * td;
            let z = if u < tau {
                Z::default()
            } else {
                let mut z = Z {
                    r: phase.cos() / (r * r),
                    i: phase.sin() / (r * r),
                };
                if receiver {
                    z = z.mul(Z::hp(f, 175e3)).mul(Z::hp(f, 350e3));
                }
                z
            };
            q.push((z.r, z.i));
        }
    }
    q
}
struct Engine {
    raw: Vec<Z>,
    valid: Vec<bool>,
    n: usize,
    rows: usize,
    per_frame: usize,
    bg: Vec<Z>,
    out: Vec<f32>,
    fs: f64,
    adc: f64,
    period: f64,
    f0: f64,
    slope: f64,
    fft_raw: Vec<f32>,
    fft_sub: Vec<f32>,
    fft_n: usize,
    event: Vec<Z>,
    noise: Vec<Vec<Z>>,
    pulses: usize,
    start: usize,
    receiver: bool,
    axes: [Vec<f64>; 3],
    scores: Vec<f32>,
    signs: Vec<f32>,
    best: usize,
    best_score: f64,
    bg_bounds: (usize, usize),
    noise_bounds: (usize, usize),
    scan_starts: Vec<usize>,
    scan_data: Vec<Vec<Z>>,
    scan_peaks: Vec<f32>,
    scan_theta: Vec<usize>,
}
impl Engine {
    fn residual(&self, k: usize) -> Vec<Z> {
        (0..self.n)
            .map(|j| {
                let z = self.raw[k * self.n + j];
                Z {
                    r: z.r - self.bg[j].r,
                    i: z.i - self.bg[j].i,
                }
            })
            .collect()
    }
    fn valid_train(&self, k: usize, l: usize) -> bool {
        k + l <= self.rows
            && k / self.per_frame == (k + l - 1) / self.per_frame
            && self.valid[k..k + l].iter().all(|v| *v)
    }
    fn train(&self, k: usize, l: usize) -> Vec<Z> {
        (k..k + l).flat_map(|a| self.residual(a)).collect()
    }
    fn spectra(&mut self, sub: bool) {
        if !(if sub { &self.fft_sub } else { &self.fft_raw }).is_empty() {
            return;
        }
        let nf = self.fft_n;
        let mut values = vec![f32::NAN; self.rows * nf];
        let weights: Vec<f64> = (0..self.n)
            .map(|j| 0.5 - 0.5 * (2. * PI * j as f64 / (self.n - 1) as f64).cos())
            .collect();
        let wsum: f64 = weights.iter().sum();
        for k in 0..self.rows {
            if !self.valid[k] {
                continue;
            }
            let mut a = vec![Z::default(); nf];
            for j in 0..self.n {
                let z = self.raw[k * self.n + j];
                let b = if sub { self.bg[j] } else { Z::default() };
                a[j] = Z {
                    r: (z.r - b.r) * weights[j],
                    i: (z.i - b.i) * weights[j],
                };
            }
            fft(&mut a);
            for j in 0..nf {
                values[k * nf + j] = (a[(j + nf / 2) % nf].power() / (wsum * wsum)) as f32;
            }
        }
        if sub {
            self.fft_sub = values
        } else {
            self.fft_raw = values
        };
    }
    fn theta(&self, index: usize) -> (f64, f64, f64) {
        let ny = self.axes[1].len();
        let nv = self.axes[2].len();
        (
            self.axes[0][index / (ny * nv)],
            self.axes[1][(index / nv) % ny],
            self.axes[2][index % nv],
        )
    }
    fn q(&self, index: usize) -> Vec<Z> {
        let (x, y, v) = self.theta(index);
        template(
            self.n,
            self.pulses,
            self.fs,
            self.adc,
            self.period,
            self.f0,
            self.slope,
            x,
            y,
            v,
            self.receiver,
        )
        .into_iter()
        .map(|(r, i)| Z { r, i })
        .collect()
    }
}
thread_local! {static ENGINE:RefCell<Option<Engine>>=const {RefCell::new(None)};}
#[no_mangle]
pub extern "C" fn allocate(n: usize) -> *mut f32 {
    let b = vec![0f32; n].into_boxed_slice();
    Box::into_raw(b) as *mut f32
}
#[no_mangle]
pub unsafe extern "C" fn release(p: *mut f32, n: usize) {
    drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(p, n)));
}
#[no_mangle]
pub unsafe extern "C" fn load(
    p: *const f32,
    rows: usize,
    n: usize,
    per_frame: usize,
    fs: f64,
    adc: f64,
    period: f64,
    f0: f64,
    slope: f64,
) {
    let data = std::slice::from_raw_parts(p, rows * n * 2);
    let raw: Vec<Z> = data
        .chunks_exact(2)
        .map(|v| Z {
            r: v[0] as f64,
            i: v[1] as f64,
        })
        .collect();
    // File has packet padding without a per-sample mask. Flag nonfinite values
    // and runs of >=8 exact complex zeros; retain isolated natural zero codes.
    let valid: Vec<bool> = raw
        .chunks_exact(n)
        .map(|row| {
            let mut run = 0;
            for z in row {
                if !z.r.is_finite() || !z.i.is_finite() {
                    return false;
                }
                if z.power() == 0. {
                    run += 1;
                    if run >= 8 {
                        return false;
                    }
                } else {
                    run = 0;
                }
            }
            true
        })
        .collect();
    ENGINE.with(|e| {
        *e.borrow_mut() = Some(Engine {
            raw,
            valid,
            n,
            rows,
            per_frame,
            bg: vec![Z::default(); n],
            out: vec![],
            fs,
            adc,
            period,
            f0,
            slope,
            fft_raw: vec![],
            fft_sub: vec![],
            fft_n: n.next_power_of_two(),
            event: vec![],
            noise: vec![],
            pulses: 1,
            start: 0,
            receiver: false,
            axes: [vec![], vec![], vec![]],
            scores: vec![],
            signs: vec![],
            best: 0,
            best_score: 0.,
            bg_bounds: (0, 0),
            noise_bounds: (0, 0),
            scan_starts: vec![],
            scan_data: vec![],
            scan_peaks: vec![],
            scan_theta: vec![],
        })
    });
}
#[no_mangle]
pub extern "C" fn set_period(period: f64) {
    ENGINE.with(|e| e.borrow_mut().as_mut().unwrap().period = period);
}
#[no_mangle]
pub extern "C" fn background(start: usize, stop: usize) -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        if start >= stop || stop > a.rows {
            return 0;
        }
        let selected: Vec<usize> = (start..stop).filter(|k| a.valid[*k]).collect();
        if selected.is_empty() {
            return 0;
        }
        a.bg_bounds = (start, stop);
        a.bg.fill(Z::default());
        for k in &selected {
            for j in 0..a.n {
                a.bg[j].r += a.raw[k * a.n + j].r / selected.len() as f64;
                a.bg[j].i += a.raw[k * a.n + j].i / selected.len() as f64;
            }
        }
        a.fft_sub.clear();
        selected.len()
    })
}
#[no_mangle]
pub extern "C" fn image(start: usize, stop: usize, mode: usize) -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        let spectral = mode >= 4;
        if spectral {
            a.spectra(mode == 5);
        }
        let height = if spectral { a.fft_n } else { a.n };
        let width = stop - start;
        a.out = vec![f32::NAN; width * height];
        for k in start..stop {
            if !a.valid[k] {
                continue;
            }
            for j in 0..height {
                let value = if spectral {
                    (if mode == 5 { &a.fft_sub } else { &a.fft_raw })[k * height + j]
                } else {
                    let z = a.raw[k * a.n + j];
                    let b = if mode == 1 || mode == 3 {
                        a.bg[j]
                    } else {
                        Z::default()
                    };
                    if mode < 2 {
                        (z.r - b.r) as f32
                    } else {
                        (z.i - b.i) as f32
                    }
                };
                a.out[j * width + k - start] = value;
            }
        }
        height
    })
}
#[no_mangle]
pub extern "C" fn trace(chirp: usize) -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        a.spectra(false);
        a.spectra(true);
        a.out = Vec::with_capacity(4 * a.n + 2 * a.fft_n);
        for mode in 0..4 {
            for j in 0..a.n {
                let z = a.raw[chirp * a.n + j];
                let b = if mode >= 2 { a.bg[j] } else { Z::default() };
                a.out.push(if mode % 2 == 0 {
                    (z.r - b.r) as f32
                } else {
                    (z.i - b.i) as f32
                });
            }
        }
        a.out
            .extend_from_slice(&a.fft_raw[chirp * a.fft_n..(chirp + 1) * a.fft_n]);
        a.out
            .extend_from_slice(&a.fft_sub[chirp * a.fft_n..(chirp + 1) * a.fft_n]);
        a.valid[chirp] as usize
    })
}
#[no_mangle]
pub extern "C" fn powers() -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        a.out = vec![f32::NAN; 3 * a.rows];
        for k in 0..a.rows {
            a.out[2 * a.rows + k] = a.valid[k] as u8 as f32;
            if a.valid[k] {
                a.out[k] = (a.raw[k * a.n..(k + 1) * a.n]
                    .iter()
                    .map(|z| z.power())
                    .sum::<f64>()
                    / a.n as f64) as f32;
                a.out[a.rows + k] =
                    (a.residual(k).iter().map(|z| z.power()).sum::<f64>() / a.n as f64) as f32;
            }
        }
        a.rows
    })
}
#[no_mangle]
pub extern "C" fn prepare(
    start: usize,
    pulses: usize,
    noise_start: usize,
    noise_stop: usize,
    receiver: usize,
) -> isize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        if ![1, 2, 4, 8, 16].contains(&pulses)
            || !a.valid_train(start, pulses)
            || noise_start >= noise_stop
            || noise_stop > a.rows
        {
            return -1;
        }
        a.start = start;
        a.noise_bounds = (noise_start, noise_stop);
        a.pulses = pulses;
        a.receiver = receiver != 0;
        a.event = a.train(start, pulses);
        let mut candidates = vec![];
        let mut k = noise_start;
        while k + pulses <= noise_stop {
            if a.valid_train(k, pulses) && (k + pulses <= start || k >= start + pulses) {
                candidates.push(k);
                k += pulses;
            } else {
                k += 1;
            }
        }
        if candidates.len() < 8 {
            return -2;
        }
        let count = candidates.len().min(24);
        a.noise = (0..count)
            .map(|i| a.train(candidates[i * candidates.len() / count], pulses))
            .collect();
        a.noise.len() as isize
    })
}
#[no_mangle]
pub extern "C" fn grid(
    xlo: f64,
    xhi: f64,
    nx: usize,
    ylo: f64,
    yhi: f64,
    ny: usize,
    vlo: f64,
    vhi: f64,
    nv: usize,
) -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        fn axis(lo: f64, hi: f64, n: usize) -> Vec<f64> {
            (0..n)
                .map(|i| {
                    if n == 1 {
                        lo
                    } else {
                        lo + (hi - lo) * i as f64 / (n - 1) as f64
                    }
                })
                .collect()
        }
        a.axes = [axis(xlo, xhi, nx), axis(ylo, yhi, ny), axis(vlo, vhi, nv)];
        a.scores = vec![f32::NAN; nx * ny * nv];
        a.signs = vec![0.; nx * ny * nv];
        a.best_score = -1.;
        a.best = 0;
        a.scan_starts.clear();
        a.scan_data.clear();
        a.scan_peaks.clear();
        a.scan_theta.clear();
        nx * ny * nv
    })
}
#[no_mangle]
pub extern "C" fn scan_range(start: usize, stop: usize, step: usize) -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        a.scan_starts = (start..stop)
            .step_by(step.max(1))
            .filter(|k| {
                a.valid_train(*k, a.pulses)
                    && (*k + a.pulses <= a.bg_bounds.0 || *k >= a.bg_bounds.1)
                    && (*k + a.pulses <= a.noise_bounds.0 || *k >= a.noise_bounds.1)
            })
            .collect();
        a.scan_data = a
            .scan_starts
            .iter()
            .map(|k| a.train(*k, a.pulses))
            .collect();
        a.scan_peaks = vec![0.; a.scan_starts.len()];
        a.scan_theta = vec![0; a.scan_starts.len()];
        a.scan_starts.len()
    })
}
#[no_mangle]
pub extern "C" fn search_batch(first: usize, count: usize) -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        let stop = (first + count).min(a.scores.len());
        for index in first..stop {
            let q = a.q(index);
            if q.iter().map(|z| z.power()).sum::<f64>() < 1e-40 {
                continue;
            }
            let mut best = 0.;
            let mut sign = 0.;
            for conj in [false, true] {
                let denom = a
                    .noise
                    .iter()
                    .map(|n| inner(&q, n, conj).power())
                    .sum::<f64>()
                    / a.noise.len() as f64;
                if denom < 1e-24 {
                    continue;
                }
                let value = inner(&q, &a.event, conj).power() / denom;
                if value > best {
                    best = value;
                    sign = conj as u8 as f32;
                }
                for (i, z) in a.scan_data.iter().enumerate() {
                    let value = (inner(&q, z, conj).power() / denom) as f32;
                    if value > a.scan_peaks[i] {
                        a.scan_peaks[i] = value;
                        a.scan_theta[i] = index;
                    }
                }
            }
            a.scores[index] = best as f32;
            a.signs[index] = sign;
            if best > a.best_score {
                a.best_score = best;
                a.best = index;
            }
        }
        stop
    })
}
#[no_mangle]
pub extern "C" fn matches() -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        a.out = a.scores.clone();
        let (x, y, v) = a.theta(a.best);
        a.out.extend_from_slice(&[
            x as f32,
            y as f32,
            v as f32,
            a.best_score as f32,
            a.signs[a.best],
            a.noise.len() as f32,
        ]);
        a.scores.len()
    })
}
#[no_mangle]
pub extern "C" fn scan_results() -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        a.out = vec![];
        for i in 0..a.scan_starts.len() {
            let (x, y, v) = a.theta(a.scan_theta[i]);
            a.out.extend_from_slice(&[
                a.scan_starts[i] as f32,
                a.scan_peaks[i],
                x as f32,
                y as f32,
                v as f32,
            ]);
        }
        a.scan_starts.len()
    })
}
#[no_mangle]
pub extern "C" fn fitted() -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        let q = a.q(a.best);
        let conj = a.signs[a.best] > 0.;
        let alpha = inner(&q, &a.event, conj);
        let norm = q.iter().map(|z| z.power()).sum::<f64>();
        a.out = vec![];
        for z in q {
            let q = if conj { Z { r: z.r, i: -z.i } } else { z };
            let fit = q.mul(Z {
                r: alpha.r / norm,
                i: alpha.i / norm,
            });
            a.out.push(fit.r as f32);
            a.out.push(fit.i as f32);
        }
        a.pulses
    })
}
#[no_mangle]
pub extern "C" fn template_values(x: f64, y: f64, v: f64, pulses: usize, receiver: usize) -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        a.out = template(
            a.n,
            pulses,
            a.fs,
            a.adc,
            a.period,
            a.f0,
            a.slope,
            x,
            y,
            v,
            receiver != 0,
        )
        .into_iter()
        .flat_map(|(r, i)| [r as f32, i as f32])
        .collect();
        a.out.len()
    })
}
#[no_mangle]
pub extern "C" fn result_ptr() -> *const f32 {
    ENGINE.with(|e| e.borrow().as_ref().unwrap().out.as_ptr())
}
#[no_mangle]
pub extern "C" fn result_len() -> usize {
    ENGINE.with(|e| e.borrow().as_ref().unwrap().out.len())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fft_signed_tone() {
        let mut a: Vec<Z> = (0..256)
            .map(|j| {
                let p = -2. * PI * 17. * j as f64 / 256.;
                Z {
                    r: p.cos(),
                    i: p.sin(),
                }
            })
            .collect();
        fft(&mut a);
        let k = a
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.power().total_cmp(&b.1.power()))
            .unwrap()
            .0;
        assert_eq!(k, 239);
        assert!((a[k].power() - 256f64.powi(2)).abs() < 1e-7);
    }
    #[test]
    fn geometry_mirror() {
        let a = template(
            225,
            4,
            12.5e6,
            2e-6,
            25.37e-6,
            77e9,
            99.9873877e12,
            0.3,
            0.1,
            341.51,
            true,
        );
        let b = template(
            225,
            4,
            12.5e6,
            2e-6,
            25.37e-6,
            77e9,
            99.9873877e12,
            -0.3,
            0.1,
            -341.51,
            true,
        );
        assert_eq!(a, b);
    }
    #[test]
    fn fitted_inner_sign() {
        let q = [Z { r: 1., i: 2. }, Z { r: 3., i: -1. }];
        let alpha = Z { r: 2., i: -3. };
        let z: Vec<Z> = q.iter().map(|q| q.mul(alpha)).collect();
        let s = inner(&q, &z, false);
        let norm = q.iter().map(|q| q.power()).sum::<f64>();
        assert!((s.r / norm - 2.).abs() < 1e-12);
        assert!((s.i / norm + 3.).abs() < 1e-12);
    }
}

#[no_mangle]
pub extern "C" fn train_values() -> usize {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        let a = e.as_mut().unwrap();
        a.out = a
            .event
            .iter()
            .flat_map(|z| [z.r as f32, z.i as f32])
            .collect();
        a.pulses
    })
}
