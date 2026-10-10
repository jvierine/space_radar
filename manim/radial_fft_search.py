"""Midpoint radial FFT search and the laboratory GUI.
Render the full scene with --disable_caching.
All plotted trajectories and phase vectors are synthetic.
"""
from pathlib import Path
import re
import h5py
import numpy as np
from manim import *
from manim_slides import Slide
from planck_to_ktb import Text
from radial_fft_assets import generate
from radial_fft_demo import generate_demo
from radial_beamforming_slides import beamforming_slides
from radial_fft_integration_slides import integration_slides
from radial_fft_cost_slides import cost_slides

FG, MUTED, BLUE, ORANGE, PURPLE = '#202830','#556575','#176BB0','#B55A00','#7939A8'

class RadialFFTSearch(Slide):
    def setup(self):
        self.camera.background_color=WHITE
        self.started=False
        self.footer=Text('juha.no/fmcw/lab/',font_size=18,color=BLUE).to_edge(DOWN,buff=.15)
        self.add(self.footer)

    def prose(self,s,size=29,color=FG):
        if re.search(r'\b(?:t0|r0|v0|a0|y0|X0|Pn|phia|lambda|gamma|theta|epsilon|fs)\b|[λφσεΔθ]',s):
            raise ValueError(f'Use eq() or texrow() for mathematical symbols: {s}')
        m=Text(s,font_size=size,color=color)
        if m.width>12.5:m.scale_to_fit_width(12.5)
        return m

    def eq(self,s,size=44):
        m=MathTex(s,font_size=size,color=FG)
        if m.width>12.5:m.scale_to_fit_width(12.5)
        return m

    def texrow(self,s,size=27):
        m=Tex(s,font_size=size,color=MUTED)
        if m.width>12.5:m.scale_to_fit_width(12.5)
        return m

    def callout(self,target,label,position,color=BLUE):
        text=self.prose(label,24,color).move_to(position)
        arrow=Arrow(text.get_bottom() if position[1]>target.get_center()[1] else text.get_top(),
                    target.get_top() if position[1]>target.get_center()[1] else target.get_bottom(),
                    color=color,buff=.12,stroke_width=3)
        self.play(FadeIn(text),GrowArrow(arrow),run_time=.6)
        self.play(Indicate(target,color=color),run_time=.6)
        return VGroup(text,arrow)

    def start(self,title,notes=''):

        # Store the completed previous slide before outgoing cleanup.
        self.next_slide(title,notes=notes)
        if self.started:
            old=[m for m in self.mobjects if m is not self.footer]
            if old:self.play(*[FadeOut(m) for m in old],run_time=.25)
        self.started=True
        self.play(FadeIn(self.prose(title,42).to_edge(UP,buff=.35)),run_time=.4)

    def content(self,equations,rows=()):
        group=VGroup(*[self.eq(s) for s in equations],*[self.prose(s,27,MUTED) if "$" not in s else self.texrow(s,27) for s in rows]).arrange(DOWN,buff=.48)
        if group.height>5.7:group.scale_to_fit_height(5.7)
        group.move_to(DOWN*.2)
        self.play(FadeIn(group),run_time=.8)
        self.wait(.7)

    def demo_image(self,name,width):
        return ImageMobject(str(Path(__file__).parent/'media/radial_fft_demo'/f'{name}.png')).scale_to_fit_width(width)

    def construct(self):
        generate_demo()
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
        times=VGroup(self.texrow(r'$R(t)$: radar-to-target distance at sample time $t$.',27),self.texrow(r'$t_0$: midpoint of the selected chirp train.',27)).arrange(DOWN,buff=.23).move_to(DOWN*2.95)
        self.play(FadeIn(times))
        self.wait(1)

        self.start('Straight-line motion still bends the range curve',
            'For constant Cartesian V, differentiate R=sqrt((Vt-x0)^2+y0^2): Rprime=V(Vt-x0)/R, Rsecond=V^2*y0^2/R^3. Valid for r0>0. These geometry identities are explanatory; the browser searches independent radial parameters rather than enforcing them.')
        radar=np.array([-5.,-.6,0]);target=np.array([-2.,1.5,0]);corner=np.array([-5.,1.5,0])
        drawing=VGroup(Dot(radar,color=BLUE),Dot(target,color=ORANGE),Line(radar,target,color=ORANGE),DashedLine(radar,corner,color=MUTED),DashedLine(corner,target,color=MUTED),Arrow([-5.8,1.5,0],[-1.1,1.5,0],buff=0,color=MUTED))
        labels=VGroup(self.prose('Radar',26,BLUE).next_to(drawing[0],DOWN),self.prose('Projectile',26,ORANGE).next_to(drawing[1],UP),self.eq(r'r_0',32).move_to([-3.25,.15,0]),self.eq(r'y_0',32).move_to([-5.5,.45,0]),self.eq(r'X_0',32).move_to([-3.6,1.1,0]))
        eqs=VGroup(self.eq(r'r_0=\sqrt{X_0^2+y_0^2}',40),self.eq(r'v_0=\frac{VX_0}{r_0}',40),self.eq(r'a_0=\frac{V^2y_0^2}{r_0^3}\geq0',40)).arrange(DOWN,buff=.5).move_to([3,.5,0])
        self.play(FadeIn(drawing),FadeIn(labels));self.play(FadeIn(eqs))
        explanation=VGroup(self.texrow(r'$X_0$: along-track separation; $y_0$: distance to the path.',26),self.texrow(r'$V$: straight-line velocity; $v_0$: radial velocity.',27),self.prose('Range can decrease while radial acceleration is positive.',27)).arrange(DOWN,buff=.4).move_to(DOWN*2.35)
        self.play(FadeIn(explanation));self.wait(1)

        self.start('Why the range, velocity and acceleration formulas are exact',
            'Constant straight-line velocity V and fixed y0. Differentiate R squared, then Rprime=VX/R. Rsecond=V squared/R minus V squared X squared/R cubed = V squared y0 squared/R cubed. Evaluate at t0. These derivatives are exact; the quadratic Taylor trajectory is an approximation away from t0.')
        self.content([
            r'X(t)=Vt-x_0,\qquad R(t)=\sqrt{X(t)^2+y_0^2}',
            r'2R\dot R=2XV\quad\Longrightarrow\quad \dot R=\frac{VX}{R}',
            r'\ddot R=\frac{V^2}{R}-\frac{V^2X^2}{R^3}=\frac{V^2y_0^2}{R^3}',
        ],[
            'A dot means differentiation with respect to time.',
            r'At $t_0$: $X=X_0$, $R=r_0$, $\dot R=v_0$, $\ddot R=a_0$.',
            r'Exact for constant $V$, fixed $y_0$ and $R>0$.',
            r'The quadratic model approximates motion away from $t_0$.',
        ])

        self.start('From a plane wave to the complex baseband signal',
            'In source-free free space Maxwell equations give the wave equation. A travelling plane wave uses exp(i omega0 t - i k0 x), with k0=omega0/c=2 pi/lambda. The physical electric field is the real part of the complex representation. Complex downconversion multiplies by exp(-i omega0 t), leaving only propagation phase. A scattering amplitude A includes polarization, path attenuation and fixed reflection phase. Subsequent round-trip animation freezes range during propagation; full FMCW phase is introduced later.')
        self.play(Write(self.eq(r'\nabla^2\mathbf E-\frac{1}{c^2}\frac{\partial^2\mathbf E}{\partial t^2}=0',43).move_to(UP*1.9)))
        plane_wave=self.eq(r'E(x,t)=A e^{i(\omega_0t-k_0x)},\qquad k_0=\frac{\omega_0}{c}=\frac{2\pi}{\lambda}',43).move_to(UP*.55)
        self.play(Write(plane_wave))
        self.play(FadeIn(self.prose('Remove the carrier oscillation by complex downconversion',27,BLUE).move_to(DOWN*.55)))
        self.play(Write(self.eq(r'z(x,t)=E(x,t)e^{-i\omega_0t}=A e^{-ik_0x}',46).move_to(DOWN*1.5)))
        self.play(FadeIn(self.texrow(r'The physical electric field is $\Re\{E\}$; $z=I+iQ$ retains amplitude and propagation phase.',26).move_to(DOWN*2.65)))
        self.wait(1)

        self.start('The echo phase comes from the round-trip distance',
            'Carrier propagation phase magnitude is wave number times total path length: (2 pi/lambda) times 2R. The outgoing and returning legs each contribute 2 pi R/lambda. Animation freezes range during propagation. With the exp(+i omega t) convention the echo has exp(-i 4 pi R/lambda); conjugated I/Q reverses this sign. A fixed reflection phase is absorbed into the complex amplitude. This is the carrier explanation; the later full FMCW model also includes chirp delay and slope.')
        self.play(Write(self.eq(r'z(t)=I(t)+iQ(t)=A e^{-ik_0\,2R(t)}=A e^{-i4\pi R(t)/\lambda}',39).move_to(UP*2.5)))
        left=np.array([-4.8,.5,0]);right=np.array([4.8,.5,0])
        endpoints=VGroup(Dot(left,color=BLUE,radius=.12),Dot(right,color=ORANGE,radius=.12))
        labels=VGroup(self.prose('Radar',28,BLUE).next_to(endpoints[0],DOWN,buff=.6),
            self.prose('Target',28,ORANGE).next_to(endpoints[1],DOWN,buff=.6))
        outgoing=Arrow(left+UP*.45,right+UP*.45,buff=.2,color=BLUE)
        returning=Arrow(right+DOWN*.45,left+DOWN*.45,buff=.2,color=ORANGE)
        self.play(FadeIn(endpoints),FadeIn(labels))
        self.play(GrowArrow(outgoing),Write(self.eq(r'\text{Outward path: }R(t)',32).move_to(UP*1.6)))
        packet=ParametricFunction(lambda u:np.array([u,.15*np.sin(15*u)*np.exp(-3*u*u),0]),
            t_range=[-1,1],color=BLUE,stroke_width=5).move_to(left+UP*.45)
        self.play(FadeIn(packet))
        self.play(packet.animate.move_to(right+UP*.45),run_time=2,rate_func=linear)
        self.play(FadeOut(packet),GrowArrow(returning))
        return_label=self.eq(r'\text{Return path: }R(t)',32).move_to(DOWN*.85)
        self.play(Write(return_label))
        packet.set_color(ORANGE).move_to(right+DOWN*.45)
        self.play(FadeIn(packet))
        self.play(packet.animate.move_to(left+DOWN*.45),run_time=2,rate_func=linear)
        self.play(FadeOut(packet))
        propagation_equation=self.eq(r'\phi_{\rm path}(t)=\underbrace{\frac{2\pi}{\lambda}}_{\text{phase per metre}}\underbrace{2R(t)}_{\text{there and back}}=\frac{4\pi R(t)}{\lambda}',36).move_to(DOWN*1.95)
        self.play(Write(propagation_equation))
        self.play(FadeIn(self.texrow(r'After carrier removal: $z(t)=A e^{-i4\pi R(t)/\lambda}=Aq(t)$, $\lambda=c/f_c$. Reflection phase is included in $A$.',23).move_to(DOWN*3.0)))
        self.wait(1)

        self.start('Coherent integration',
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
        self.play(FadeIn(self.texrow(r'All echo samples now have the same phase: their complex sum is $NA$.',27).move_to(DOWN*2.95)))
        self.play(FadeIn(self.texrow(r'$q$: predicted echo; $A$: amplitude; $N$: acquired samples. Noiseless illustration.',22).move_to(DOWN*3.4)))
        self.wait(1)


        self.start('Least squares exposes the matched-filter inner product',
            'For complex z=Aq(theta)+n, expand J=(z-Aq)^H(z-Aq). The cross term contains q^H z, the model-measurement inner product and temporal matched filter.')
        self.play(Write(self.eq(r'z=Aq(\theta)+n,\qquad J(A,\theta)=\|z-Aq(\theta)\|^2',41).move_to(UP*1.8)))
        expansion=MathTex(r'J=',r'z^Hz',r'+|A|^2q^Hq',r'-2\operatorname{Re}\{A^*',r'q^Hz',r'\}',font_size=43,color=FG).move_to(UP*.3)
        self.play(Write(expansion))
        self.callout(expansion[1],'Measured energy',[-4.7,1.3,0],MUTED)
        self.callout(expansion[2],'Predicted energy',[-1,-1.05,0],ORANGE)
        self.callout(expansion[4],'Model–measurement inner product',[4,-1.05,0],BLUE)
        self.play(Write(self.eq(r's(\theta)=q(\theta)^Hz=\sum_jq_j(\theta)^*z_j',44).set_color(BLUE).move_to(DOWN*2)))
        self.play(FadeIn(self.texrow(r'$z$: complex measurements; $q$: predicted echo; $A$: complex amplitude; $n$: noise.',25).move_to(DOWN*2.95)))
        self.play(FadeIn(self.texrow(r'$H$: conjugate transpose; $j$: sample over the entire acquired train; $\theta=(r_0,v_0,a_0)$.',24).move_to(DOWN*3.4)))
        self.wait(1)

        self.start('Fit the amplitude; maximize normalized matched energy',
            'Complete the square with Q=q^Hq>0 and s=q^Hz. Minimizing over A gives Ahat=s/Q and Jmin=z^Hz-|s|^2/Q. Therefore trajectory fitting maximizes normalized matched energy. White Gaussian noise makes least squares likelihood fitting; independent temporal noise gives matched variance Pn Q. Correct-template SNR is |A|^2 Q/Pn, N times per-sample SNR for unit-magnitude samples. Colored noise requires inverse covariance weighting.')
        self.play(Write(self.eq(r'Q=q^Hq,\quad s=q^Hz,\quad J=Q\left|A-\frac{s}{Q}\right|^2+z^Hz-\frac{|s|^2}{Q}',39).move_to(UP*1.9)))
        self.play(Write(self.eq(r'\widehat A=\frac{s}{Q},\qquad \widehat\theta=\operatorname*{arg\,max}_{\theta}\frac{|q(\theta)^Hz|^2}{q(\theta)^Hq(\theta)}',41).move_to(UP*.6)))
        self.play(Write(self.eq(r'\mathbb E|q^Hn|^2=P_nQ,\qquad \rho=\frac{|s|^2}{P_nQ}',38).move_to(DOWN*.55)))
        self.play(Write(self.eq(r'\mathrm{SNR}_{\rm signal,out}=\frac{|A|^2Q}{P_n}=N\frac{|A|^2}{P_n}\quad\text{if }|q_j|=1',36).move_to(DOWN*1.65)))
        self.play(FadeIn(self.texrow(r'$Q$: template energy; $P_n$: noise power per sample; $N$: acquired samples.',25).move_to(DOWN*2.65)))
        self.play(FadeIn(self.texrow(r'With colored noise: $s=q^HC_n^{-1}z$, $Q=q^HC_n^{-1}q$; $C_n$: temporal noise covariance.',24).move_to(DOWN*3.2)))
        self.wait(1)

        self.start('Acceleration adds a changing phase',
            'Synthetic uniform carrier-phase example at 77 GHz. Isolates the acceleration term for teaching; production uses full FMCW phase, including fast-time chirp terms. Convention q=exp(-i4pi R/lambda).')
        self.play(Write(VGroup(self.eq(r'h=t-t_0,\quad \Delta R_a=\tfrac12a_0h^2,\quad \phi_a=-\frac{2\pi a_0h^2}{\lambda}',38),self.eq(r'c_j=e^{i\phi_a(h_j)}\quad\text{acceleration factor at sample }j',32)).arrange(DOWN,buff=.35).move_to(UP*2.3)))
        self.play(FadeIn(self.eq(r'h:\ \text{time from midpoint};\quad\lambda:\ \text{wavelength};\quad\phi_a:\ \text{acceleration phase}',26).move_to(UP*1.15)))
        self.play(FadeIn(self.demo_image('acceleration_factor',12.5).move_to(DOWN*.95)))
        self.wait(1)

        self.start('Undo the trial acceleration in the complex voltage',
            'Multiply measured complex voltage by conjugate of the trial acceleration factor. Correct trial removes chirping; wrong trials leave changing phase. A unit-magnitude factor changes phase without changing instantaneous power.')
        self.play(FadeIn(self.eq(r'z_{\rm corrected,j}=z_jc_j^*\quad\text{(undo the acceleration phase)}',38).move_to(UP*2.6)))
        self.play(FadeIn(self.demo_image('complex_voltage',12).move_to(DOWN*.2)))
        self.play(FadeIn(self.prose('Conjugation reverses phase. The correct trial leaves a constant-frequency wave.',25,MUTED).move_to(DOWN*3.05)))
        self.wait(1)

        self.start('One FFT scores all velocity bins for this correction',
            'Uniform-sampling teaching example: Doppler frequency f=-2v/lambda. FFT evaluates all sampled-frequency inner products at once. In full FMCW, fast frequency mixes range and velocity and slow frequency carries Doppler; repeat bounded correction groups rather than claiming one FFT covers every physical template.')
        self.play(FadeIn(self.demo_image('velocity_search',12).move_to(UP*.1)))
        self.play(FadeIn(self.prose('Try acceleration → correct the voltage → FFT → read every velocity-bin power.',27).move_to(DOWN*2.45)))
        self.play(FadeIn(self.eq(r'f_D=-\frac{2v}{\lambda}\quad\text{where }v\text{ is radial velocity}',35).move_to(DOWN*3.05)))
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
        self.play(FadeIn(self.texrow(r'$N$: acquired samples; $N_v$: velocity bins; $L$: padded FFT length.',23).move_to(DOWN*3.4)))
        self.wait(1)

        self.start('An FFT is the matched-filter sum, evaluated for every bin',
            'Exact algebra for factored sampled templates at DFT bins: q_mj=c_j exp(i2pi mj/L), s_m=sum z_j conjugate(q_mj)=DFT of z conjugate(c) at bin m. N acquired samples, L FFT length with zero padding; m bin, j sample. Exact factorization is not proof that approximate grouped production corrections match every exact physical template; direct verification follows.')
        fft_math = VGroup(*[self.eq(line,36) for line in [
            r'\text{From the acceleration slide: }c_j=e^{i\phi_a(h_j)},\quad\phi_a(h_j)=-\frac{2\pi a_0h_j^2}{\lambda}',
            r'q_{m,j}=\underbrace{c_j}_{\text{acceleration}}\underbrace{e^{2\pi i m j/L}}_{\text{velocity-bin sinusoid}}',
            r's_m=\sum_{j=0}^{N-1}z_jq_{m,j}^*=\sum_{j=0}^{N-1}(z_jc_j^*)e^{-2\pi i m j/L}',
            r's_m=\operatorname{FFT}(z\,c^*)_m',
        ]], *[self.eq(line,27) for line in [
            r'z_j:\ \text{measured voltage};\quad c_j^*:\ \text{undoes the earlier acceleration phase}',
            r'q_{m,j}:\ \text{predicted echo};\quad j:\ \text{sample};\quad m:\ \text{frequency bin}',
            r'N:\ \text{acquired samples};\quad L:\ \text{padded FFT length}',
        ]], self.prose('Zero-pad before the FFT: each bin is one matched-filter sum.',25,MUTED)).arrange(DOWN,buff=.22)
        if fft_math.height>5.7:fft_math.scale_to_fit_height(5.7)
        fft_math.move_to(DOWN*.1)
        self.play(FadeIn(fft_math));self.wait(1)


        self.start('Coarse search: complex FFT pseudocode',
            'For both possible stored complex phase conventions, loop over the automatic correction groups. Each group contains one acceleration and a bounded velocity interval. Correct complex voltages, perform fast FFT per chirp, then slow FFT on the physical grid column union. Preserve complex outputs. Map every physical node through signed, modulo-wrapped frequency bins, score four receiver energies by expected matched noise. Keep distinct alias branches and full-bank peaks. Final GLS fit is on the following slide.')
        code=[
            (0,r'\textbf{grid} $\gets\operatorname{Grid}(\text{bounds},L_{\rm dB})$',FG),
            (0,r'\textbf{for} $\sigma\in\{+1,-1\}$: \quad \textit{complex phase convention}',MUTED),
            (1,r'\textbf{for each} correction group $g$:',FG),
            (2,r'$c_{g,k,j}\gets e^{i\sigma\psi_g(k,j)}$',ORANGE),
            (2,r'$\mathcal M_g\gets\operatorname{UniqueFastBins}(\text{grid}_g,\sigma)$',FG),
            (2,r'\textbf{for each} receiver $b$:',FG),
            (3,r'$x_{b,k,j}\gets z_{b,k,j}c_{g,k,j}^{*}$',ORANGE),
            (3,r'$F_{b,k,:}\gets\operatorname{FFT}_{j}(x_{b,k,:})$',BLUE),
            (3,r'$S_{b,:,m}\gets\operatorname{FFT}_{k}(F_{b,:,m}),\quad m\in\mathcal M_g$',BLUE),
            (2,r'\textbf{for each} physical node $\theta=(r_0,v_0,a_0)\in\text{grid}_g$:',FG),
            (3,r'$(m,\ell)\gets\operatorname{WrappedBins}(\theta,\sigma)$',FG),
            (3,r'$\rho(\theta,\sigma)\gets\sum_b|S_{b,\ell,m}|^2/(Q\sum_bP_{n,b})$',FG),
            (0,r'\textbf{seeds} $\gets$ bank peaks and strongest candidates per alias branch',PURPLE),
        ]
        rows=VGroup()
        for indent,line,color in code:
            row=self.texrow(line,25).set_color(color)
            row.to_edge(LEFT,buff=.65+indent*.36)
            rows.add(row)
        rows.arrange(DOWN,buff=.18,aligned_edge=LEFT)
        # Preserve indentation after vertically arranging the rows.
        for row,(indent,_,_) in zip(rows,code):row.to_edge(LEFT,buff=.65+indent*.36)
        rows.move_to(DOWN*.05)
        self.play(LaggedStart(*[FadeIn(row) for row in rows],lag_ratio=.12))
        self.play(FadeIn(self.texrow(r'$Q=\sum_{k,j}|q_{k,j}|^2$; $g$: correction group; $m,\ell$: fast and slow bins.',22).move_to(DOWN*3.35)))
        self.wait(1)

        self.start('Final fit: noise-weighted refinement and coherent beamforming',
            'Exact full-FMCW templates score the shortlisted physical seeds by profiled joint GLS. Starting from the strongest seed, bounded Nelder-Mead optimizes the three motion parameters while solving four complex amplitudes analytically at every trial. Per-sample receive covariance C comes from quiet residuals, temporally white assumption. The same fitted amplitudes give w proportional C^-1 Ahat and complex beam sBF=w^H s. Normalize weights without changing SNR; propagate Q w^H C w. Alternative aliases are retained, not declared uniquely resolved.')
        code=[
            r'\textbf{function} $\operatorname{FitScore}(\theta,\sigma)$:',
            r'\qquad $q\gets\operatorname{FullFMCWEcho}(\theta,\sigma)$',
            r'\qquad $\mathbf s\gets\sum_jq_j^*\mathbf z_j,\quad Q\gets\sum_j|q_j|^2$',
            r'\qquad $\widehat{\mathbf A}\gets\mathbf s/Q$',
            r'\qquad \textbf{return} $\mathbf s^H C^{-1}\mathbf s/Q$',
            r'$(\theta_0,\sigma_*)\gets\operatorname*{arg\,max}_{\rm seeds}\operatorname{FitScore}(\theta,\sigma)$',
            r'$\theta_*\gets\operatorname{BoundedNelderMeadMax}(\operatorname{FitScore},\theta_0,\sigma_*)$',
            r'$\mathbf s_*,Q_*,\widehat{\mathbf A}_*\gets$ full-template evaluation at $\theta_*$',
            r'$\mathbf w\gets C^{-1}\widehat{\mathbf A}_*/\|C^{-1}\widehat{\mathbf A}_*\|$',
            r'$s_{\rm BF}\gets\mathbf w^H\mathbf s_*,\quad P_{n,\rm BF}\gets Q_*\mathbf w^H C\mathbf w$',
        ]
        rows=VGroup(*[self.texrow(line,27).set_color(PURPLE if i>=5 else FG) for i,line in enumerate(code)]).arrange(DOWN,buff=.25,aligned_edge=LEFT)
        if rows.width>12:rows.scale_to_fit_width(12)
        rows.move_to(DOWN*.1)
        self.play(LaggedStart(*[FadeIn(row) for row in rows],lag_ratio=.15))
        self.play(FadeIn(self.texrow(r'$C$: noise covariance; $\widehat{\mathbf A}$: fitted RX amplitudes; $\mathbf w$: beam weights.',22).move_to(DOWN*3.25)))
        self.wait(1)

        self.start('The acceleration correction in the full FMCW echo',
            'Exact phase factorization at a representative correction group: phiFMCW=2pi[-(f0+gamma u) tau + gamma tau squared/2], tau=2R/c. R=r0+v0h+a0h squared/2. Residual psi subtracts dominant fast/slow sinusoids, including chirp resets and delay squared term. Phase-only c=exp(i psi), multiplied conjugate before FFT. Carrier-only acceleration contribution reduces to -2pi a0h squared/lambda from the earlier slide. Shared corrections are approximate away from representative parameters; final templates use full phase.')
        self.content([
            r'R=r_0+v_0h+\tfrac12a_0h^2,\quad\tau=\frac{2R}{c},\quad h=t-t_0',
            r'\phi_{\rm FMCW}=2\pi\left[-(f_0+\gamma u)\tau+\tfrac12\gamma\tau^2\right]',
            r'\psi=\phi_{\rm FMCW}-2\pi[f_{\rm fast}(u-\bar u)+f_{\rm slow}d]',
            r'c_{k,j}=e^{i\psi_{k,j}},\quad x_{b,k,j}=z_{b,k,j}c_{k,j}^*',
            r'\text{Earlier acceleration contribution: }\phi_a=-\frac{2\pi a_0h^2}{\lambda}',
        ],[
            'Remove residual phase; leave the two sinusoids for the FFTs to search.',
            r'$u$: within-chirp time; $\bar u$: midpoint; $d$: chirp-start offset.',
            r'$f_0$: start frequency; $\gamma$: slope; $c$: speed of light.',
            r'$\tau$: delay; $\lambda$: wavelength; $b$: RX; $k$: chirp; $j$: sample.',
        ])

        integration_slides(self)

        self.start('You give the bounds. The software builds the grid.',
            'Default bounds r0 .001..3 m, v0 0..600 m/s, a0 0..1e6 m/s2. FFT bin spacing sets r/v resolution; global phase derivative bounds set acceleration spacing and velocity groups. The full modeled coherent integration loss includes physical grid spacing, shared corrections and FFT-bin lookup. Oversized grids are rejected without truncating bounds; 64 million-node cap.')
        bounds=VGroup(self.eq(r'0.001\leq r_0\leq3\ \mathrm{m}',44),self.eq(r'0\leq v_0\leq600\ \mathrm{m/s}',44),self.eq(r'0\leq a_0\leq10^6\ \mathrm{m/s^2}',44)).arrange(DOWN,buff=.6).move_to(UP*.3)
        self.play(FadeIn(bounds))
        self.play(FadeIn(self.prose('Grid: trial combinations of range, velocity and acceleration.',29).move_to(DOWN*2.1)))
        self.play(FadeIn(self.prose('Allowed loss: tolerated reduction of matched signal power.',28,MUTED).move_to(DOWN*2.9)))
        self.wait(1)

        self.start('Choose the loss of the full coherent integration',
            'The normalized full-template overlap is the retained matched signal power. For unit magnitude samples its loss is at most the population variance of the phase mismatch, after removing constant phase. The Rust grid budgets the RMS contributions of physical parameter-grid spacing, correction sharing and rounded FFT-bin lookup over every acquired sample. Default 2 dB, not an endpoint phase bound. Carrier-only acceleration illustrates why acceleration grid count grows as elapsed train duration squared; production uses full FMCW derivative bounds.')
        overlap=MathTex(r'\eta=',r'\frac{|\langle q_{\rm true},q_{\rm bank}\rangle|^2}{\|q_{\rm true}\|^2\|q_{\rm bank}\|^2}',font_size=43,color=FG).move_to(UP*1.65)
        self.play(Write(overlap))
        self.callout(overlap[1],'Retained power after the entire coherent sum',[0,2.8,0],BLUE)
        self.play(Write(self.eq(r'\ell=1-\eta\leq\operatorname{Var}_j(\delta\phi_j),\qquad\ell_{\max}=1-10^{-L_{\rm dB}/10}',37).move_to(UP*.25)))
        self.play(Write(self.eq(r'L_{\rm dB}=2\quad\Longrightarrow\quad\eta\geq10^{-2/10}\approx0.631',38).set_color(BLUE).move_to(DOWN*.65)))
        self.play(Write(self.eq(r'\delta\phi_j=-\frac{2\pi\delta a}{\lambda}h_j^2,\quad\sigma_{\phi}=\frac{2\pi|\delta a|}{\lambda}\operatorname{SD}_j(h_j^2)',35).move_to(DOWN*1.6)))
        self.play(Write(self.eq(r'|\delta a|\leq\Delta a/2,\quad N_a-1=\left\lceil\frac{a_{\max}-a_{\min}}{\Delta a}\right\rceil\propto T_{\rm span}^2',35).move_to(DOWN*2.45)))
        self.play(FadeIn(self.texrow(r'$h_j=t_j-t_0$; $\operatorname{SD}_j$: spread over all acquired samples; $T_{\rm span}$: elapsed train duration.',22).move_to(DOWN*3.1)))
        self.wait(1)

        cost_slides(self)

        self.start('Estimate noise from the quiet beat signal',
            'Per-RX residual noise variance uses every intact selected quiet chirp and all fast-time samples. Quiet mean is independently fitted per RX and fast-time sample. Division by Ns*(M-1) is the sample-variance correction; repeated stationary wall power is removed. One quiet chirp cannot separate arbitrary repeated echoes from noise: code uses a raw-power upper bound and omits RCS.')
        image=ImageMobject(str(Path(__file__).parent/'assets/quiet_noise_example.png')).scale_to_fit_width(11.7).move_to(UP*1.75)
        self.play(FadeIn(image))
        self.play(FadeIn(self.prose('Use the blue Background window, before decoding or FFT.',25,BLUE).move_to(ORIGIN)))
        formula=MathTex(r'P_{n,b}=',r'\frac{\sum_{k,j}|z^{\rm bg}_{b,k,j}|^2}{N_s(M-1)}',font_size=42,color=FG).move_to(DOWN*1.1)
        self.play(Write(formula))
        noise_arrow=self.callout(formula[1],'Squared magnitude of the residual voltage',[-3.8,-2.55,0],ORANGE)
        noise_arrow[1].put_start_and_end_on(noise_arrow[1].get_start(),formula[1].get_corner(UL)+RIGHT*.2+DOWN*.1)
        self.play(FadeIn(self.texrow(r'$k$: quiet chirp; $j$: sample; $b$: receiver.',25).move_to([3.3,-2.5,0])))
        self.play(FadeIn(self.texrow(r'$N_s$: samples per chirp; $M$: intact background chirps.',25).move_to(DOWN*3.15)))
        self.play(FadeIn(self.texrow(r'$M-1$ corrects for the mean fitted from these same chirps.',22).move_to(DOWN*3.5)))
        self.wait(1)

        self.start('The matched filter also reduces noise',
            'Assume temporally white residual noise. Unnormalized matched sum has noise energy Pn*sum|q|^2, which is the denominator in the next slide. For a phase-only template, |q|=1 and averaging N samples has variance Pn/N. Full complex sample-rate bandwidth fs scales to analysis bandwidth 1/Tcoh. Acquired sample time excludes idle/frame gaps.')
        self.play(FadeIn(self.prose('For the matched sum, expected noise power is',30).move_to(UP*2.1)))
        self.play(Write(self.eq(r'P_{\rm noise,sum}=P_{n,b}\sum_j|q_j|^2',46).move_to(UP*1.2)))
        self.play(FadeIn(self.texrow(r'For a phase-only template averaged over $N$ samples,',30).move_to(UP*.25)))
        self.play(Write(self.eq(r'P_{\rm noise,average}=\frac{P_{n,b}}{N}=\frac{P_{n,b}}{f_sT_{\rm coh}}',42).move_to(DOWN*.45)))
        terms=VGroup(self.texrow(r'$N$: acquired samples used in coherent integration.',28),self.eq(r'f_s:\ \text{sample rate};\quad T_{\rm coh}=\frac{N}{f_s}:\ \text{acquired coherent time}',32),self.eq(r'\text{Noise bandwidth: }f_s\longrightarrow\frac{1}{T_{\rm coh}}',32)).arrange(DOWN,buff=.32).move_to(DOWN*2.45)
        self.play(FadeIn(terms));self.wait(1)

        self.start('Then add the four matched powers',
            'rho is total observed matched power / expected matched noise power, not excess-signal SNR. The noise denominator assumes temporal white noise, independently estimated per RX after subtracting stationary quiet mean. Displayed SNR subtracts expected noise and floors its dB display at zero; physical powers retain raw values.')
        self.play(Write(self.eq(r'\rho=\frac{\sum_b|s_b|^2}{(\sum_bP_{n,b})(\sum_j|q_j|^2)}',54).move_to(UP*1.8)))
        description=VGroup(self.prose('Top: add the matched powers of the four receivers.',30,BLUE),self.prose('Bottom: expected noise power after the same filtering.',30,ORANGE)).arrange(DOWN,buff=.35).move_to(UP*.1)
        self.play(FadeIn(description))
        terms=VGroup(self.eq(r'|s_b|^2:\ \text{matched power from receiver }b',32),self.eq(r'P_{n,b}:\ \text{background noise power per sample, receiver }b',32),self.eq(r'\sum_j|q_j|^2:\ \text{energy of the predicted template}',32),self.eq(r'\rho:\ \text{matched power divided by expected noise power}',32)).arrange(DOWN,buff=.35).move_to(DOWN*1.95)
        self.play(FadeIn(terms));self.wait(1)

        self.start('Search the 3D array; project it only for display',
            'The discrete search maximum is argmax over the full three-dimensional range velocity acceleration bank. The RV and VA maps are only maximum projections for visualization and do not define separate searches. Candidate checks and continuous GLS refine the discrete search result.')
        self.play(Write(self.eq(r'(r_*,v_*,a_*)=\operatorname*{arg\,max}_{r_0,v_0,a_0}\rho(r_0,v_0,a_0)',43).move_to(UP*2)))
        layers=VGroup()
        for k in range(3):
            layer=VGroup(*[Rectangle(width=.42,height=.42,stroke_color=MUTED,fill_color=BLUE,fill_opacity=.15).move_to([-4.8+i*.48+k*.15,.9-j*.48+k*.15,0]) for j in range(3) for i in range(4)])
            layers.add(layer)
        self.play(LaggedStart(*[FadeIn(l) for l in layers],lag_ratio=.3))
        peak=layers[1][6];self.play(peak.animate.set_fill(ORANGE,opacity=1),Indicate(peak,color=ORANGE))
        self.play(FadeIn(self.prose('Largest voxel in the full 3D bank',23,ORANGE).move_to([-3.9,-1.1,0])))
        maps=VGroup(self.eq(r'M_{rv}(r_0,v_0)=\max_{a_0}\rho(r_0,v_0,a_0)',36),
            self.eq(r'M_{va}(v_0,a_0)=\max_{r_0}\rho(r_0,v_0,a_0)',36)).arrange(DOWN,buff=.6).move_to([2.8,.1,0])
        self.play(GrowArrow(Arrow([-1.6,.2,0],[-.5,.2,0],color=BLUE)),FadeIn(maps))
        self.play(FadeIn(self.prose('Two 2D pictures of that same array—not two new searches.',27,BLUE).move_to(DOWN*2.15)))
        self.play(FadeIn(self.prose('The selected 3D candidate starts the joint motion / receiver-amplitude fit.',25).move_to(DOWN*2.85)))
        self.wait(1)

        beamforming_slides(self)

        self.start('From matched SNR to received signal power',
            'Thermal calibration assumption Tsys=9000 K; white noise across complex sample-rate bandwidth. rho is observed power/noise, so signal SNR S=max(rho-1,0). B=1/Tcoh with acquired sample time. Do not invert the GUI zero-dB display floor. Joint beam SNR removes fitted-amplitude noise bias nu, and received power divides by covariance-calibrated combining gain.')
        self.content([
            r'S_{\rm RX}=\max(\rho-1,0),\quad S_{\rm beam}=\max(\rho_{\rm GLS}-\nu,0),\quad B=\frac{1}{T_{\rm coh}}',
            r'P_r=S k_B T_{\rm sys} B=\frac{S k_B T_{\rm sys}}{T_{\rm coh}}',
        ],[
            r'$S$: linear signal-to-noise power ratio, after subtracting expected noise.',
            r'$B$: noise bandwidth; $T_{\mathrm{coh}}$: acquired integration time.',
            r'$k_B$: Boltzmann constant; $T_{\mathrm{sys}}$: assumed noise temperature, $9000\ \mathrm{K}$.',
            r'$P_r$: received power in watts; use the underlying ratio before the display floor.',
        ])

        self.start('Use distance to convert received power to radar cross section',
            'Monostatic radar equation with linear gains and total linear loss L. Assumes far-field monostatic geometry and stated transmitter/receiver calibration. Effective R is fitted range or explicit override. Individual RX results use own noise and matched power; beam uses the implemented noise-weighted combining gain with equal calibrated antenna-response assumptions.')
        self.content([
            r'\sigma=\frac{P_r(4\pi)^3R^4 L}{P_tG_tG_r\lambda^2}',
        ],[
            r'$\sigma$: radar cross section in $\mathrm{m^2}$; $R$: radar-to-projectile distance.',
            r'$P_t$: transmitted power; $G_t,G_r$: linear transmit and receive gains.',
            r'$\lambda$: wavelength; $L$: total loss factor ($1$ means no loss).',
            'Compute separately for each receiver, using its own SNR and noise estimate.',
            'For the beam, divide by the fitted covariance-calibrated combining gain.',
            'Absolute RCS depends on the assumed temperature, gains, losses and distance.',
        ])

        self.start('Find every metallic-sphere diameter consistent with that RCS',
            'Mie theory gives the electromagnetic scattering of an ideal perfectly conducting sphere. Compare predicted monostatic cross section with inferred cross section; all crossings inside the user diameter bounds are reported. Multiple diameters can produce one cross section. Same inversion per RX and beam result. Synthetic 77GHz illustration; not a measured projectile.')
        self.play(FadeIn(self.demo_image('mie_diameter',11.8).move_to(UP*.25)))
        self.play(FadeIn(self.prose('Mie theory predicts RCS versus diameter; each crossing is a possible diameter.',27).move_to(DOWN*2.55)))
        self.play(FadeIn(self.texrow(r'Candidate envelope: $0.949$–$2.277\ \mathrm{mm}$; keep all three roots.',26).move_to(DOWN*3.1)))
        self.wait(1)

        self.start('Select background. Select analysis. Press Play.',
            'Real GUI screenshot supplied by Juha Vierinen. Blue Background selects quiet chirps used for mean and residual noise estimation. Purple Analyze selects the train scan interval. Yellow coherent integration marker is the current acquired chirp train. Press Play to analyze successive intact trains. Drag windows to move, edges to resize; empty-space drag zooms; Reset full view restores the viewport. Download HDF5 exports completed estimates and current bank.')
        image=ImageMobject(str(Path(__file__).parent/'assets/lab_gui_example.png')).scale_to_fit_width(11.1).move_to(DOWN*.2)
        self.play(FadeIn(image))
        def region(x,y,width,height,color):
            return Rectangle(width=image.width*width,height=image.height*height,stroke_color=color,stroke_width=5).move_to(image.get_corner(DL)+RIGHT*image.width*x+UP*image.height*y)
        background=region(.492,.32,.082,.395,BLUE)
        analysis=region(.66,.32,.235,.395,PURPLE)
        play=region(.66,.915,.06,.072,ORANGE)
        caption=self.prose('1. Blue Background: quiet chirps for mean removal and noise.',25,BLUE).move_to(DOWN*3.15)
        self.play(Create(background),FadeIn(caption));self.wait(.7)
        self.play(Create(analysis),Transform(caption,self.prose('2. Purple Analyze: choose the interval containing the target.',25,PURPLE).move_to(DOWN*3.15)));self.wait(.7)
        self.play(Create(play),Transform(caption,self.prose('3. Press Play. Yellow marks the chirp train being integrated.',25,ORANGE).move_to(DOWN*3.15)))
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
