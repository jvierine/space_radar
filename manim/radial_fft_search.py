"""Midpoint radial FFT search and the laboratory GUI.
Render the full scene with --disable_caching.
All plotted trajectories and phase vectors are synthetic.
"""
import h5py
import numpy as np
from manim import *
from manim_slides import Slide
from planck_to_ktb import Text
from radial_fft_assets import generate

FG, MUTED, BLUE, ORANGE, PURPLE = '#202830','#556575','#176BB0','#B55A00','#7939A8'

class RadialFFTSearch(Slide):
    def setup(self):
        self.camera.background_color=WHITE
        self.started=False
        self.footer=Text('juha.no/fmcw/lab/',font_size=18,color=BLUE).to_edge(DOWN,buff=.15)
        self.add(self.footer)

    def prose(self,s,size=29,color=FG):
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

    def construct(self):
        path=generate()
        with h5py.File(path) as f:psi=f['correction_rad'][:,112]

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

        self.start('Align the chirps before adding them',
            'Synthetic eight-chirp phasors after removing the dominant fast/slow sinusoid. Each direction is remaining phase. Correcting that phase restores coherent addition. Unsampled gaps retain phase timing. This is for multi-chirp integration.')
        circles=VGroup(*[Circle(radius=1.25,color=MUTED).move_to([x,.35,0]) for x in [-3,3]])
        arrows=VGroup(*[Arrow([-3,.35,0],[-3+1.15*np.cos(p),.35+1.15*np.sin(p),0],buff=0,color=ORANGE,stroke_width=3) for p in psi])
        aligned=VGroup(*[Arrow([3,.35,0],[4.15,.35,0],buff=0,color=BLUE,stroke_width=3) for _ in psi])
        labels=VGroup(self.prose('Different phases',32,ORANGE).move_to([-3,-1.35,0]),self.prose('Aligned phases',32,BLUE).move_to([3,-1.35,0]))
        self.play(Create(circles),FadeIn(arrows),FadeIn(labels))
        self.play(TransformFromCopy(arrows,aligned),run_time=1.2)
        self.play(FadeIn(self.prose('Coherent integration: align phases, then add the signals.',32).move_to(DOWN*2.35)))
        self.play(FadeIn(VGroup(self.prose('A chirp is one transmitted frequency sweep.',26,MUTED),self.prose('Each arrow shows phase: the angle of a received chirp signal.',26,MUTED)).arrange(DOWN,buff=.2).move_to(DOWN*3)))
        self.wait(1)

        self.start('Try a correction, then let the FFT search',
            'Actual algorithm: acceleration planes and bounded velocity correction groups, representative range-dependent residual phase, fast FFT for each RX/chirp then sparse slow FFT. Search explicit r0/v0/a0 nodes and both I/Q orientations; incoherently add the four RX matched powers. Verify candidates with exact quadratic templates and refine. Fast frequency mixes range and Doppler; correction contains all parameters.')
        factor=self.eq(r'\text{signal}=\text{dominant sinusoid}\times\text{phase correction}',41).move_to(UP*1.7)
        self.play(FadeIn(factor))
        boxes=VGroup()
        for x,label,color in [(-4.3,'Try a correction',ORANGE),(0,'FFT search',BLUE),(4.3,'Keep the best fit',PURPLE)]:
            box=RoundedRectangle(width=3.65,height=1.15,corner_radius=.12,color=color)
            box.add(self.prose(label,30,color)).move_to([x,0,0]);boxes.add(box)
        connectors=VGroup(*[Arrow([x,0,0],[x+.6,0,0],buff=0,color=MUTED) for x in [-2.35,1.75]])
        self.play(LaggedStart(*[FadeIn(b) for b in boxes],lag_ratio=.3),FadeIn(connectors))
        rows=VGroup(self.prose('Sinusoid: a constant-frequency oscillation.',28),self.prose('Correction: removes extra phase changes predicted by the model.',28),self.prose('FFT: fast Fourier transform; searches frequencies efficiently.',28),self.prose('Repeat corrections and keep the strongest match.',28)).arrange(DOWN,buff=.3).move_to(DOWN*2)
        self.play(FadeIn(rows));self.wait(1)

        self.start('You give the bounds. The software builds the grid.',
            'Default bounds r0 .001..3 m, v0 0..600 m/s, a0 0..1e6 m/s2. FFT bin spacing sets r/v resolution; global phase derivative bounds set acceleration spacing and velocity groups. The phase tolerance is stagewise, not a bound on total model error. Oversized grids are rejected without truncating bounds; 32 million-node cap.')
        bounds=VGroup(self.eq(r'0.001\leq r_0\leq3\ \mathrm{m}',44),self.eq(r'0\leq v_0\leq600\ \mathrm{m/s}',44),self.eq(r'0\leq a_0\leq10^6\ \mathrm{m/s^2}',44)).arrange(DOWN,buff=.6).move_to(UP*.3)
        self.play(FadeIn(bounds))
        self.play(FadeIn(self.prose('Grid: trial combinations of range, velocity and acceleration.',29).move_to(DOWN*2.1)))
        self.play(FadeIn(self.prose('Allowed loss: tolerated reduction of matched signal power.',28,MUTED).move_to(DOWN*2.9)))
        self.wait(1)

        self.start('First, match each receiver to the predicted signal',
            'The inner product is s_b=sum_j conjugate(q_j)*z_bj over acquired samples. It is evaluated separately per RX and for both stored-I/Q orientations. The template can include the optional receiver response at direct verification.')
        self.play(Write(self.eq(r's_b=\sum_j q_j^*z_{b,j}',58).move_to(UP*1.7)))
        terms=VGroup(
            self.eq(r'q_j:\ \text{predicted complex signal at sample }j',34),
            self.eq(r'z_{b,j}:\ \text{background-subtracted data from receiver }b',34),
            self.eq(r'q_j^*:\ \text{complex conjugate; reverses the predicted phase}',34),
            self.eq(r'\sum_j:\ \text{add all samples in the selected chirp train}',34),
            self.eq(r's_b:\ \text{the resulting complex sum for receiver }b',34),
        ).arrange(DOWN,buff=.38).move_to(DOWN*.7)
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
