"""Midpoint radial FFT search and the laboratory GUI.
Render the full scene with --disable_caching.
All plotted trajectories and phase vectors are synthetic.
"""
from pathlib import Path
import h5py
import numpy as np
from manim import *
from manim_slides import Slide
from planck_to_ktb import Text
from radial_fft_assets import generate
from radial_fft_demo import generate_demo
from radial_beamforming_slides import beamforming_slides
from radial_fft_integration_slides import integration_slides

FG, MUTED, BLUE, ORANGE, PURPLE = '#202830','#556575','#176BB0','#B55A00','#7939A8'

class RadialFFTSearch(Slide):
    def setup(self):
        self.camera.background_color=WHITE
        self.started=False
        self.footer=Text('juha.no/fmcw/lab/',font_size=18,color=BLUE).to_edge(DOWN,buff=.15)
        self.add(self.footer)

    def prose(self,s,size=29,color=FG):
        for symbol,word in [('λ','lambda'),('φ','phi'),('σ','sigma'),('ε','epsilon'),('Δ','Delta '),('⁻¹²','^(-12)'),('i²','i squared')]:s=s.replace(symbol,word)
        m=Text(s,font_size=size,color=color)
        if m.width>12.5:m.scale_to_fit_width(12.5)
        return m

    def eq(self,s,size=44):
        m=MathTex(s,font_size=size,color=FG)
        if m.width>12.5:m.scale_to_fit_width(12.5)
        return m

    def start(self,title,notes=''):
        # Store the completed previous slide before outgoing cleanup.
        self.next_slide(title,notes=notes)
        if self.started:
            old=[m for m in self.mobjects if m is not self.footer]
            if old:self.play(*[FadeOut(m) for m in old],run_time=.25)
        self.started=True
        self.play(FadeIn(self.prose(title,42).to_edge(UP,buff=.35)),run_time=.4)

    def content(self,equations,rows=()):
        group=VGroup(*[self.eq(s) for s in equations],*[self.prose(s,27,MUTED) for s in rows]).arrange(DOWN,buff=.48)
        if group.height>5.7:group.scale_to_fit_height(5.7)
        group.move_to(DOWN*.2)
        self.play(FadeIn(group),run_time=.8)
        self.wait(.7)

    def demo_image(self,name,width):
        return ImageMobject(str(Path(__file__).parent/'media/radial_fft_demo'/f'{name}.png')).scale_to_fit_width(width)

    def construct(self):
        demo=generate_demo()
        path=generate()
        with h5py.File(path) as f:
            phase=f['phase_rad'][:];beat_times=(f['h_s'][:]-f['h_s'][:].min())*1e6

        self.start('Coherent integration and projectile diameter',
            'Two goals: maximize signal-to-noise ratio by matching motion across chirps; infer an ideal metallic sphere diameter from calibrated received power. The laboratory uses a full FMCW phase model.')
        cover=VGroup(
            self.prose('How can we maximize signal-to-noise ratio',37,BLUE),
            self.prose('with coherent integration?',37,BLUE),
            self.prose('How can we estimate projectile diameter',37,ORANGE),
            self.prose('from the measured signal-to-noise ratio?',37,ORANGE),
        ).arrange(DOWN,buff=.4).move_to(UP*.1)
        self.play(FadeIn(cover))
        self.play(FadeIn(self.prose('Motion model → phase correction → FFT matched filter → diameter',27,MUTED).move_to(DOWN*2.7)))
        self.wait(1)

        self.start('Describe the motion with three numbers',
            'Parameters are defined at the acquired train midpoint t0. v0 is signed radial velocity, not Cartesian speed; a0 is radial range curvature. The quadratic range model is local.')
        model=self.eq(r'R(t)=r_0+v_0(t-t_0)+\tfrac12a_0(t-t_0)^2',54).move_to(UP*1.5)
        definitions=VGroup(
            self.eq(r'r_0:\ \text{distance at the midpoint}',39),
            self.eq(r'v_0:\ \text{how fast the distance changes}',39),
            self.eq(r'a_0:\ \text{how that velocity changes}',39),
        ).arrange(DOWN,buff=.55).move_to(DOWN*.7)
        self.play(Write(model));self.play(FadeIn(definitions))
        times=VGroup(self.prose('R(t) is radar-to-target distance at sample time t.',27,MUTED),self.prose('t0 is the midpoint of the selected chirp train.',27,MUTED)).arrange(DOWN,buff=.23).move_to(DOWN*2.95)
        self.play(FadeIn(times))
        self.wait(1)

        self.start('Straight-line motion still bends the range curve',
            'For constant Cartesian V, differentiate R=sqrt((Vt-x0)^2+y0^2): Rprime=V(Vt-x0)/R, Rsecond=V^2*y0^2/R^3. Valid for r0>0. These geometry identities are explanatory; the browser searches independent radial parameters rather than enforcing them.')
        radar=np.array([-5.,-.6,0]);target=np.array([-2.,1.5,0]);corner=np.array([-5.,1.5,0])
        drawing=VGroup(Dot(radar,color=BLUE),Dot(target,color=ORANGE),Line(radar,target,color=ORANGE),DashedLine(radar,corner,color=MUTED),DashedLine(corner,target,color=MUTED),Arrow([-5.8,1.5,0],[-1.1,1.5,0],buff=0,color=MUTED))
        labels=VGroup(self.prose('Radar',26,BLUE).next_to(drawing[0],DOWN),self.prose('Projectile',26,ORANGE).next_to(drawing[1],UP),self.eq(r'r_0',32).move_to([-3.25,.15,0]),self.eq(r'y_0',32).move_to([-5.5,.45,0]),self.eq(r'X_0',32).move_to([-3.6,1.1,0]))
        eqs=VGroup(self.eq(r'r_0=\sqrt{X_0^2+y_0^2}',40),self.eq(r'v_0=\frac{VX_0}{r_0}',40),self.eq(r'a_0=\frac{V^2y_0^2}{r_0^3}\geq0',40)).arrange(DOWN,buff=.5).move_to([3,.5,0])
        self.play(FadeIn(drawing),FadeIn(labels));self.play(FadeIn(eqs))
        explanation=VGroup(self.prose('X0: along-track separation. y0: distance from radar to the path.',26),self.prose('V: straight-line velocity. v0: radial velocity.',27),self.prose('Range can decrease while radial acceleration is positive.',27)).arrange(DOWN,buff=.4).move_to(DOWN*2.35)
        self.play(FadeIn(explanation));self.wait(1)

        self.start('Why the range, velocity and acceleration formulas are exact',
            'Constant straight-line velocity V and fixed y0. Differentiate R squared, then Rprime=VX/R. Rsecond=V squared/R minus V squared X squared/R cubed = V squared y0 squared/R cubed. Evaluate at t0. These derivatives are exact; the quadratic Taylor trajectory is an approximation away from t0.')
        self.content([
            r'X(t)=Vt-x_0,\qquad R(t)=\sqrt{X(t)^2+y_0^2}',
            r'2R\dot R=2XV\quad\Longrightarrow\quad \dot R=\frac{VX}{R}',
            r'\ddot R=\frac{V^2}{R}-\frac{V^2X^2}{R^3}=\frac{V^2y_0^2}{R^3}',
        ],[
            'A dot means differentiation with respect to time.',
            'At t0: X = X0, R = r0, velocity = v0, acceleration = a0.',
            'Exact for constant V, fixed y0 and nonzero range.',
            'The quadratic model approximates motion away from t0.',
        ])

        self.start('Write down the coherent integration',
            'The inner product is s_b=sum_j conjugate(q_j)*z_bj over acquired samples. It is evaluated separately per RX and for both stored-I/Q orientations. The template can include the optional receiver response at direct verification.')
        self.play(FadeIn(self.prose('To add the echo constructively: undo its predicted phase, then sum.',28).move_to(UP*2.35)))
        self.play(Write(self.eq(r's_b=\sum_j q_j^*z_{b,j}',46).move_to(UP*1.3)))
        terms=VGroup(
            self.eq(r'q_j:\ \text{predicted complex signal at sample }j',30),
            self.eq(r'z_{b,j}:\ \text{background-subtracted data from receiver }b',30),
            self.eq(r'q_j^*:\ \text{complex conjugate; reverses the predicted phase}',30),
            self.eq(r'\sum_j:\ \text{add all samples in the selected chirp train}',30),
            self.eq(r's_b:\ \text{the resulting complex sum for receiver }b',30),
        ).arrange(DOWN,buff=.38).move_to(DOWN*1.25)
        self.play(FadeIn(terms));self.wait(1)

        self.start('Unwind the echo waveform, then add it constructively',
            'Synthetic full FMCW beat phase: transmitted chirp phase evaluated at the delayed return minus local transmit phase, with delay 2R(t)/c and quadratic moving range. Each acquired sample is multiplied by conjugate of its full predicted phase. Animation continuously removes that phase; phase unwrapping alone would not change the voltage or align it. No data in idle gaps. Noise would remain random rather than flatten.')
        self.play(FadeIn(self.prose('The beat phase contains the transmit sweep, propagation delay and target motion.',26,MUTED).move_to(UP*2.5)))
        axes=Axes(x_range=[0,float(beat_times.max()),25],y_range=[-1,1,.5],x_length=12,y_length=2.9,
            axis_config={'color':MUTED,'include_tip':False}).move_to(UP*.2)
        self.play(Create(axes))
        self.play(FadeIn(self.prose('Time across eight chirps (µs); gaps contain no samples',25,MUTED).move_to(DOWN*1.65)))
        p=phase-phase[0,0]
        def curves(fraction):
            z=np.exp(1j*(1-fraction)*p)
            group=VGroup()
            for k in range(8):
                for values,col in [(z[k].real,BLUE),(z[k].imag,ORANGE)]:
                    line=VMobject(color=col,stroke_width=2)
                    line.set_points_as_corners([axes.c2p(x,y) for x,y in zip(beat_times[k],values)])
                    group.add(line)
            return group
        wave=curves(0)
        self.play(FadeIn(wave))
        legend=VGroup(self.prose('Real voltage',25,BLUE),self.prose('Imaginary voltage',25,ORANGE)).arrange(RIGHT,buff=.7).move_to(UP*1.95)
        self.play(FadeIn(legend));self.wait(1)
        correction=self.eq(r'z_jq_j^*=A\quad\text{when }z_j=Aq_j',40).move_to(DOWN*2.25)
        self.play(FadeIn(correction))
        self.play(UpdateFromAlphaFunc(wave,lambda m,alpha:m.become(curves(alpha))),run_time=5,rate_func=linear)
        self.play(FadeIn(self.prose('All echo samples now have the same phase: their complex sum is N × A.',27).move_to(DOWN*2.95)))
        self.play(FadeIn(self.prose('q: full predicted echo. A: echo amplitude. N: acquired samples. Synthetic noiseless signal.',22,MUTED).move_to(DOWN*3.4)))
        self.wait(1)

        self.start('Why the matched filter maximizes signal-to-noise ratio',
            'For z=Aq+n and white complex noise variance Pn per acquired sample, output w dagger z has signal power |A| squared |w dagger q| squared and variance Pn ||w|| squared. Cauchy-Schwarz gives upper bound |A| squared ||q|| squared/Pn, reached for w proportional q. For phase-only templates ||q|| squared=N. Colored noise requires inverse covariance weighting; the GUI assumes temporal white noise after quiet subtraction.')
        self.content([
            r'z=Aq+n,\qquad \mathrm{SNR}_{\rm out}=\frac{|A|^2|w^Hq|^2}{P_n\|w\|^2}',
            r'|w^Hq|^2\leq\|w\|^2\|q\|^2\quad\Longrightarrow\quad w\propto q',
            r'\mathrm{SNR}_{\rm max}=N\frac{|A|^2}{P_n}\quad\text{for }|q_j|=1',
        ],[
            'z: acquired voltage vector. q: predicted echo. n: noise. A: echo amplitude.',
            'w: integration weights. H: conjugate transpose. Squared norm: sum of powers.',
            'Pn: noise variance per sample. N: acquired samples.',
            'The inequality is Cauchy–Schwarz; matching the echo attains its bound.',
            'Assumes independent, equal-variance noise across time samples.',
        ])

        self.start('Acceleration adds a changing phase',
            'Synthetic uniform carrier-phase example at 77 GHz. Isolates the acceleration term for teaching; production uses full FMCW phase, including fast-time chirp terms. Convention q=exp(-i4pi R/lambda).')
        self.play(Write(self.eq(r'h=t-t_0,\quad \Delta R_a=\tfrac12a_0h^2,\quad \phi_a=-\frac{2\pi a_0h^2}{\lambda}',40).move_to(UP*2)))
        self.play(FadeIn(self.prose('h: time from midpoint. λ: wavelength. φa: acceleration phase.',23,MUTED).move_to(UP*1.25)))
        self.play(FadeIn(self.demo_image('acceleration_factor',12.5).move_to(DOWN*.75)))
        self.play(FadeIn(self.prose('Complex factor = cos(φa) + i sin(φa); i² = -1.',27).move_to(DOWN*2.8)))
        self.play(FadeIn(self.prose('Synthetic Doppler example; the GUI uses the full FMCW phase.',23,MUTED).move_to(DOWN*3.25)))
        self.wait(1)

        self.start('Undo the trial acceleration in the complex voltage',
            'Multiply measured complex voltage by conjugate of the trial acceleration factor. Correct trial removes chirping; wrong trials leave changing phase. A unit-magnitude factor changes phase without changing instantaneous power.')
        self.play(FadeIn(self.prose('Corrected voltage = measured voltage × conjugate trial factor',30).move_to(UP*2.6)))
        self.play(FadeIn(self.demo_image('complex_voltage',12).move_to(DOWN*.2)))
        self.play(FadeIn(self.prose('Conjugation reverses phase. The correct trial leaves a constant-frequency wave.',25,MUTED).move_to(DOWN*3.05)))
        self.wait(1)

        self.start('One FFT scores all velocity bins for this correction',
            'Uniform-sampling teaching example: Doppler frequency f=-2v/lambda. FFT evaluates all sampled-frequency inner products at once. In full FMCW, fast frequency mixes range and velocity and slow frequency carries Doppler; repeat bounded correction groups rather than claiming one FFT covers every physical template.')
        self.play(FadeIn(self.demo_image('velocity_search',12).move_to(UP*.1)))
        self.play(FadeIn(self.prose('Try acceleration → correct the voltage → FFT → read every velocity-bin power.',27).move_to(DOWN*2.45)))
        self.play(FadeIn(self.prose('Here Doppler frequency f = -2v / λ; v is radial velocity.',26,MUTED).move_to(DOWN*3.05)))
        self.wait(1)

        self.start('How can we do this efficiently?',
            'For each acceleration correction perform an FFT, rather than one direct N-sample inner product per frequency. An L-point FFT costs O(L log2 L), including zero padding. Direct evaluation at Nv velocity bins costs O(N Nv). N is acquired samples; L is transform length. The diagram isolates acceleration/Doppler; full FMCW adds bounded correction groups and range mapping.')
        self.play(FadeIn(self.prose('Loop over acceleration corrections. Each FFT scores many velocities together.',29).move_to(UP*2.5)))
        branches=VGroup()
        for row,(label,col) in enumerate([(r'a=0',MUTED),(r'a=\Delta a',ORANGE),(r'a=2\Delta a',PURPLE)]):
            y=1.35-row*1.15
            trial=self.eq(label,29).set_color(col).move_to([-5.4,y,0])
            corr=RoundedRectangle(width=3.2,height=.7,corner_radius=.08,color=col).move_to([-2.7,y,0])
            corr.add(self.prose('Conjugate correction',24,col).move_to(corr))
            fft=RoundedRectangle(width=1.35,height=.7,corner_radius=.08,color=BLUE).move_to([.25,y,0])
            fft.add(self.prose('FFT',28,BLUE).move_to(fft))
            arrow=Arrow(corr.get_right(),fft.get_left(),buff=.12,color=MUTED)
            bars=VGroup(*[Rectangle(width=.1,height=.12+.48*np.exp(-((k-8-row)/2)**2),fill_color=BLUE,fill_opacity=.7,stroke_width=0).align_to([0,y-.3,0],DOWN) for k in range(20)]).arrange(RIGHT,buff=.035,aligned_edge=DOWN).move_to([3.5,y,0])
            branches.add(VGroup(trial,corr,arrow,fft,Arrow(fft.get_right(),bars.get_left(),buff=.12,color=MUTED),bars))
        self.play(LaggedStart(*[FadeIn(row) for row in branches],lag_ratio=.35))
        self.play(FadeIn(self.prose('Each bar is a velocity-bin matched power; keep the best across all rows.',26).move_to(DOWN*2.2)))
        self.play(Write(self.eq(r'\text{Direct: }O(NN_v)\qquad\text{FFT: }O(L\log_2L)',35).move_to(DOWN*2.85)))
        self.play(FadeIn(self.prose('N: acquired samples. Nv: velocity bins. L: FFT length (including padding).',23,MUTED).move_to(DOWN*3.4)))
        self.wait(1)

        self.start('An FFT is the matched-filter sum, evaluated for every bin',
            'Exact algebra for factored sampled templates at DFT bins: q_mj=c_j exp(i2pi mj/L), s_m=sum z_j conjugate(q_mj)=DFT of z conjugate(c) at bin m. N acquired samples, L FFT length with zero padding; m bin, j sample. Exact factorization is not proof that approximate grouped production corrections match every exact physical template; direct verification follows.')
        self.content([
            r'q_{m,j}=c_j e^{2\pi i m j/L}',
            r's_m=\sum_{j=0}^{N-1}z_jq_{m,j}^*=\sum_{j=0}^{N-1}(z_jc_j^*)e^{-2\pi i m j/L}',
            r's_m=\operatorname{FFT}(z\,c^*)_m',
        ],[
            'zj: complex measurement. cj: trial phase correction. q: predicted signal.',
            'j: sample index. m: frequency bin. N: samples. L: FFT length.',
            'Zero-pad to L samples; each FFT bin is one matched-filter sum.',
            'Two FFT axes apply the same identity within and between chirps.',
        ])

        self.start('The same answers, with shared computation',
            'Independent NumPy direct phasor sums and FFT agree at all 356 velocity bins, for all three acceleration trials. Timings are medians of 101 runs, one BLAS thread, precomputed direct phasors; common correction excluded. This is a synthetic CPU illustration, not a browser benchmark.')
        self.play(FadeIn(self.demo_image('fft_equals_direct',8).move_to([-2,.15,0])))
        bench=VGroup(
            self.prose(f"{demo['velocity_bins']} velocity bins",29,BLUE),
            self.prose(f"Direct: {demo['cpu_direct_ms']:.3f} ms",27),
            self.prose(f"FFT: {demo['cpu_fft_ms']:.3f} ms",27),
            self.prose(f"{demo['speedup']:.1f}× faster in this example",26,BLUE),
            self.prose('All bin powers agree',26),
            self.prose('to within 1e-12.',26),
        ).arrange(DOWN,buff=.35).move_to([4.2,.2,0])
        self.play(FadeIn(bench))
        self.play(FadeIn(self.prose('Synthetic NumPy CPU comparison; precomputed direct templates, 101-run median.',23,MUTED).move_to(DOWN*3.05)))
        self.wait(1)

        self.start('What our Rust / WebGPU implementation actually does',
            'Source: lab-core/src/fft_search.rs transform/search_cell, web/lab/gpu.mjs. Each acceleration plane has bounded velocity correction groups. Fast FFT per chirp, sparse slow FFT over selected frequency columns; physical node mapping handles aliases within user bounds. Add four RX powers for each I/Q orientation. Verify candidate nodes with full quadratic FMCW templates and refine. CPU fallback.')
        rows=VGroup(*[self.prose(line,29,col) for line,col in [
            ('1. Build correction groups across the allowed acceleration and velocity bounds.',FG),
            ('2. For each receiver: multiply samples by the conjugate correction.',ORANGE),
            ('3. FFT within each chirp; then FFT across chirps at needed frequency bins.',BLUE),
            ('4. Map FFT bins to range / velocity / acceleration; add four RX powers.',FG),
            ('5. Normalize by matched noise; keep candidates and refine with full templates.',PURPLE),
        ]]).arrange(DOWN,buff=.55).move_to(UP*.1)
        self.play(LaggedStart(*[FadeIn(row) for row in rows],lag_ratio=.3))
        self.play(FadeIn(self.prose('Repeat for both I/Q orientations. WebGPU when available; Rust/Wasm fallback.',24,MUTED).move_to(DOWN*2.8)))
        self.play(FadeIn(self.prose('Correction sharing is approximate; final candidates use the full FMCW model.',24,MUTED).move_to(DOWN*3.25)))
        self.wait(1)

        integration_slides(self)

        self.start('You give the bounds. The software builds the grid.',
            'Default bounds r0 .001..3 m, v0 0..600 m/s, a0 0..1e6 m/s2. FFT bin spacing sets r/v resolution; global phase derivative bounds set acceleration spacing and velocity groups. The phase tolerance is stagewise, not a bound on total model error. Oversized grids are rejected without truncating bounds; 32 million-node cap.')
        bounds=VGroup(self.eq(r'0.001\leq r_0\leq3\ \mathrm{m}',44),self.eq(r'0\leq v_0\leq600\ \mathrm{m/s}',44),self.eq(r'0\leq a_0\leq10^6\ \mathrm{m/s^2}',44)).arrange(DOWN,buff=.6).move_to(UP*.3)
        self.play(FadeIn(bounds))
        self.play(FadeIn(self.prose('Grid: trial combinations of range, velocity and acceleration.',29).move_to(DOWN*2.1)))
        self.play(FadeIn(self.prose('Allowed loss: tolerated reduction of matched signal power.',28,MUTED).move_to(DOWN*2.9)))
        self.wait(1)

        self.start('How many acceleration corrections are needed?',
            'Carrier-phase bound with monostatic radial acceleration and midpoint origin: adjacent trial acceleration spacing Delta a gives maximum phase separation 2pi Delta a H squared/lambda where H=max_j|t_j-t0|. Require separation <=epsilon. Na=1+ceil((amax-amin)/Delta a) includes both endpoints. At nearest grid point acceleration mismatch <=Delta a/2. Fixed bounds and epsilon give Na-1 proportional H squared. In hard_target gmf_opts.py the corresponding start-referenced total-path convention uses pi Delta A tau squared/lambda with A=2a. GUI uses full FMCW phase-derivative bounds, not only carrier formula.')
        ax=Axes(x_range=[-1,1,.5],y_range=[0,4,1],x_length=4.4,y_length=2.4,axis_config={'color':MUTED,'include_tip':False}).move_to([-4.3,.9,0])
        curves=VGroup(*[ax.plot(lambda x,k=k:k*x*x,x_range=[-1,1],color=col) for k,col in enumerate([BLUE,ORANGE,PURPLE,MUTED])])
        self.play(Create(ax),Create(curves))
        self.play(FadeIn(self.prose('Phase separation / epsilon',22,MUTED).move_to([-4.3,2.5,0])))
        self.play(FadeIn(self.prose('Adjacent phase curves stay within ε',22).move_to([-4.3,-.7,0])))
        self.play(FadeIn(self.prose('Farthest samples: h = -H and +H',22,MUTED).move_to([-4.3,-1.2,0])))
        formulas=VGroup(
            self.eq(r'|\Delta\phi|_{\max}=\frac{2\pi\Delta a\,H^2}{\lambda}\leq\epsilon',34),
            self.eq(r'\Delta a\leq\frac{\epsilon\lambda}{2\pi H^2}',34),
            self.eq(r'N_a=1+\left\lceil\frac{a_{\max}-a_{\min}}{\Delta a}\right\rceil',34),
        ).arrange(DOWN,buff=.55).move_to([2.5,.4,0])
        self.play(FadeIn(formulas))
        terms=VGroup(self.prose('Δa: grid spacing. ε: allowed adjacent phase separation (radians).',24),self.prose('H: largest sample time offset. λ: wavelength. Na: trial count (round up).',24),self.prose('Double the elapsed train span → about four times as many accelerations.',25,BLUE),self.prose('For gapped chirps use elapsed span here; noise T(coh) counts acquired sample time.',22,MUTED)).arrange(DOWN,buff=.22).move_to(DOWN*2.65)
        self.play(FadeIn(terms));self.wait(1)

        self.start('Estimate noise from the quiet beat signal',
            'Per-RX residual noise variance uses every intact selected quiet chirp and all fast-time samples. Quiet mean is independently fitted per RX and fast-time sample. Division by Ns*(M-1) is the sample-variance correction; repeated stationary wall power is removed. One quiet chirp cannot separate arbitrary repeated echoes from noise: code uses a raw-power upper bound and omits RCS.')
        steps=VGroup(self.prose('For each receiver, subtract its quiet mean waveform.',31),self.prose('Use all selected quiet samples: real squared + imaginary squared.',29)).arrange(DOWN,buff=.35).move_to(UP*1.65)
        self.play(FadeIn(steps))
        self.play(Write(self.eq(r'P_{n,b}=\frac{\sum_{k,j}|z^{\rm bg}_{b,k,j}|^2}{N_s(M-1)}',52).move_to(UP*.1)))
        terms=VGroup(self.prose('z(bg): background-subtracted complex beat signal in the quiet window.',27),self.prose('b: receiver. k: quiet chirp. j: sample within a chirp.',27),self.prose('Ns: samples per chirp. M: intact quiet chirps (at least two).',27),self.prose('M - 1 corrects for estimating the quiet mean from those chirps.',27)).arrange(DOWN,buff=.27).move_to(DOWN*2.05)
        self.play(FadeIn(terms));self.wait(1)

        self.start('The matched filter also reduces noise',
            'Assume temporally white residual noise. Unnormalized matched sum has noise energy Pn*sum|q|^2, which is the denominator in the next slide. For a phase-only template, |q|=1 and averaging N samples has variance Pn/N. Full complex sample-rate bandwidth fs scales to analysis bandwidth 1/Tcoh. Acquired sample time excludes idle/frame gaps.')
        self.play(FadeIn(self.prose('For the matched sum, expected noise power is',30).move_to(UP*2.1)))
        self.play(Write(self.eq(r'P_{\rm noise,sum}=P_{n,b}\sum_j|q_j|^2',46).move_to(UP*1.2)))
        self.play(FadeIn(self.prose('For a phase-only template averaged over N samples,',30).move_to(UP*.25)))
        self.play(Write(self.eq(r'P_{\rm noise,average}=\frac{P_{n,b}}{N}=\frac{P_{n,b}}{f_sT_{\rm coh}}',46).move_to(DOWN*.7)))
        terms=VGroup(self.prose('N: acquired samples used in the coherent integration.',28),self.prose('fs: sample rate. T(coh) = N / fs: acquired coherent time.',28),self.prose('Noise bandwidth changes from fs to 1 / T(coh).',28)).arrange(DOWN,buff=.3).move_to(DOWN*2.3)
        self.play(FadeIn(terms));self.wait(1)

        self.start('Then add the four matched powers',
            'rho is total observed matched power / expected matched noise power, not excess-signal SNR. The noise denominator assumes temporal white noise, independently estimated per RX after subtracting stationary quiet mean. Displayed SNR subtracts expected noise and floors its dB display at zero; physical powers retain raw values.')
        self.play(Write(self.eq(r'\rho=\frac{\sum_b|s_b|^2}{(\sum_bP_{n,b})(\sum_j|q_j|^2)}',54).move_to(UP*1.8)))
        description=VGroup(self.prose('Top: add the matched powers of the four receivers.',30,BLUE),self.prose('Bottom: expected noise power after the same filtering.',30,ORANGE)).arrange(DOWN,buff=.35).move_to(UP*.1)
        self.play(FadeIn(description))
        terms=VGroup(self.eq(r'|s_b|^2:\ \text{matched power from receiver }b',32),self.eq(r'P_{n,b}:\ \text{background noise power per sample, receiver }b',32),self.eq(r'\sum_j|q_j|^2:\ \text{energy of the predicted template}',32),self.eq(r'\rho:\ \text{matched power divided by expected noise power}',32)).arrange(DOWN,buff=.35).move_to(DOWN*1.95)
        self.play(FadeIn(terms));self.wait(1)

        self.start('Each map pixel keeps the strongest match',
            'Mrv is MAX over acceleration at fixed r0/v0; Mva is MAX over range at fixed v0/a0. rho is the three-dimensional bank score. The map color is its monotonic display-SNR transform; ties can result from the zero dB floor.')
        first=VGroup(self.eq(r'M_{rv}(r_0,v_0)=\max_{a_0}\rho(r_0,v_0,a_0)',43),self.prose('Range versus velocity: keep the best acceleration at each pixel.',29)).arrange(DOWN,buff=.4).move_to(UP*1.15)
        second=VGroup(self.eq(r'M_{va}(v_0,a_0)=\max_{r_0}\rho(r_0,v_0,a_0)',43),self.prose('Velocity versus acceleration: keep the best range at each pixel.',29)).arrange(DOWN,buff=.4).move_to(DOWN*1.05)
        self.play(FadeIn(first));self.play(FadeIn(second))
        self.play(FadeIn(self.prose('M is the map score. MAX selects the largest score.',29,MUTED).move_to(DOWN*2.8)))
        self.wait(1)

        beamforming_slides(self)

        self.start('From matched SNR to received signal power',
            'Thermal calibration assumption Tsys=9000 K; white noise across complex sample-rate bandwidth. rho is observed power/noise, so signal SNR S=max(rho-1,0). B=1/Tcoh with acquired sample time. Do not invert the GUI zero-dB display floor. Beam received-power estimate additionally divides by assumed ideal four-RX gain.')
        self.content([
            r'S=\max(\rho-1,0),\qquad B=\frac{1}{T_{\rm coh}}',
            r'P_r=S k_B T_{\rm sys} B=\frac{S k_B T_{\rm sys}}{T_{\rm coh}}',
        ],[
            'S: signal-to-noise power ratio in linear units, after subtracting expected noise.',
            'B: analysis noise bandwidth. T(coh): acquired integration time.',
            'kB: Boltzmann constant. T(sys): assumed system noise temperature, 9000 K.',
            'Pr: received signal power in watts. Use the underlying ratio, not the display floor.',
        ])

        self.start('Use distance to convert received power to radar cross section',
            'Monostatic radar equation with linear gains and total linear loss L. Assumes far-field monostatic geometry and stated transmitter/receiver calibration. Effective R is fitted range or explicit override. Individual RX results use own noise and matched power; beam uses ideal four-channel gain, approximate with unequal/correlated RX.')
        self.content([
            r'\sigma=\frac{P_r(4\pi)^3R^4 L}{P_tG_tG_r\lambda^2}',
        ],[
            'σ: radar cross section, in square metres. R: radar-to-projectile distance.',
            'Pt: transmitted power. Gt and Gr: transmit and receive gains, in linear units.',
            'λ: radar wavelength. L: total loss factor (1 means no loss).',
            'Compute separately for each receiver, using its own SNR and noise estimate.',
            'For four-RX beamforming, divide out the assumed ideal gain of four.',
            'Absolute RCS depends on the assumed temperature, gains, losses and distance.',
        ])

        self.start('Find every metallic-sphere diameter consistent with that RCS',
            'Mie theory gives the electromagnetic scattering of an ideal perfectly conducting sphere. Compare predicted monostatic cross section with inferred cross section; all crossings inside the user diameter bounds are reported. Multiple diameters can produce one cross section. Same inversion per RX and beam result. Synthetic 77GHz illustration; not a measured projectile.')
        self.play(FadeIn(self.demo_image('mie_diameter',11.8).move_to(UP*.25)))
        self.play(FadeIn(self.prose('Mie theory predicts RCS versus diameter; each crossing is a possible diameter.',27).move_to(DOWN*2.55)))
        self.play(FadeIn(self.prose('Use diameter bounds or other evidence to choose among multiple solutions.',26,MUTED).move_to(DOWN*3.1)))
        self.wait(1)

        self.start('Select background. Select analysis. Press Play.',
            'GUI schematic, not a screenshot. Blue and purple windows move and resize. Drag empty space to zoom; Reset full view restores the recording. Moving yellow coherent train analyzes that train. Quiet mean and residual noise estimated independently per RX from all intact selected background samples. B=1/Tcoh, Tcoh=N/fs. WebGPU has CPU fallback.')
        panel=Rectangle(width=12,height=2.5,color=MUTED).shift(UP*.5)
        bg=Rectangle(width=2,height=2.5,color=BLUE,fill_color=BLUE,fill_opacity=.14).move_to([-3.7,.5,0])
        scan=Rectangle(width=5.4,height=2.5,color=PURPLE,fill_color=PURPLE,fill_opacity=.10).move_to([1.4,.5,0])
        coh=Rectangle(width=.45,height=2.5,color=ORANGE,fill_color=ORANGE,fill_opacity=.2).move_to([-.7,.5,0])
        self.play(FadeIn(panel))
        self.play(FadeIn(bg),FadeIn(self.prose('1. Background',30,BLUE).move_to([-3.7,2.3,0])))
        self.play(FadeIn(scan),FadeIn(coh),FadeIn(self.prose('2. Analysis',30,PURPLE).move_to([1.4,2.3,0])))
        play=RoundedRectangle(width=2.2,height=.7,corner_radius=.1,color=BLUE).move_to([0,-1.5,0]);play.add(self.prose('3. Play',32,BLUE).move_to(play))
        self.play(FadeIn(play))
        self.play(FadeIn(self.prose('Move or resize the coloured windows to choose the intervals.',29).move_to(DOWN*2.6)))
        self.wait(1)

        self.start('Watch the results appear',
            'Maps are MAX projections v0 versus r0 and v0 versus a0. Waveforms/continuous fitted range are lines, scan histories are scatter. Receiver angles are relative phases, not arrival directions. RCS/diameters depend on calibration assumptions, sphere inversion can yield multiple roots. A fit alone does not prove detection. Completed HDF5 includes scalar history, phase outputs and current bank/fit; URL carries GUI state only.')
        cards=VGroup()
        for x,title,body,color in [(-4.3,'Motion',['Range','Radial velocity','Radial acceleration'],BLUE),(0,'Receivers',['Phase: signal angle','SNR: signal-to-noise ratio'],PURPLE),(4.3,'Target',['RCS: radar cross section','Effective scattering area','Sphere diameter'],ORANGE)]:
            card=RoundedRectangle(width=3.8,height=2,corner_radius=.12,color=color).move_to([x,.5,0])
            card.add(self.prose(title,36,color).move_to([x,1.05,0]))
            item=VGroup(*[self.prose(line,23) for line in body]).arrange(DOWN,buff=.16)
            if item.width>3.4:item.scale_to_fit_width(3.4)
            card.add(item.move_to([x,.05,0]));cards.add(card)
        self.play(LaggedStart(*[FadeIn(card) for card in cards],lag_ratio=.3))
        rows=VGroup(self.prose('Results update as the analysis window is scanned.',31),self.prose('HDF5 (Hierarchical Data Format): all results in one file.',31),self.prose('Share the URL to share your current view.',31)).arrange(DOWN,buff=.5).move_to(DOWN*2)
        self.play(FadeIn(rows));self.wait(1)
        # Preserve the completed final slide; no outgoing cleanup.
