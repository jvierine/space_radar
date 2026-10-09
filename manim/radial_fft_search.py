"""Midpoint radial FFT search and the laboratory GUI.
Render the full scene with --disable_caching. SHOW_PROVENANCE=0 hides footers.
All plotted trajectories and phase vectors are synthetic.
"""
import os
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
        self.footer=Text('Sources: radial_fft_search.py; radial_fft_assets.py; lab-core/src/radial_search.rs',font_size=13,color=MUTED).to_edge(DOWN,buff=.12)
        if os.getenv('SHOW_PROVENANCE','1')!='0':self.add(self.footer)

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
        self.start('Range, radial velocity and range curvature')
        eq=self.eq(r'R(t)\simeq r_0+v_0h+\tfrac12a_0h^2,\qquad h=t-t_\star').move_to(UP*1.9)
        ax=Axes(x_range=[-100,100,50],y_range=[.76,.84,.02],x_length=9,y_length=2.8,axis_config={'color':MUTED,'include_numbers':True,'font_size':22,'decimal_number_config':{'color':FG}},tips=False).shift(DOWN*.25)
        curve=ax.plot(lambda x:.8-300*x*1e-6+.5*5e5*(x*1e-6)**2,color=ORANGE)
        labels=VGroup(self.prose('Time from midpoint (microseconds)',23).next_to(ax,DOWN),self.prose('Slant range (m)',23).next_to(ax,LEFT).rotate(PI/2))
        note=self.prose('Synthetic example: negative radial velocity, positive range curvature',24,MUTED).to_edge(DOWN,buff=.65)
        self.play(FadeIn(eq),Create(ax),Create(curve),FadeIn(labels),FadeIn(note));self.wait(.7)

        self.start('Fast time, chirp time and the midpoint','Chirp period is reconstructed; gaps enter phase timing, not acquired signal energy.')
        bars=VGroup(*[Rectangle(width=1.10,height=.6,color=BLUE,fill_color=BLUE,fill_opacity=.18).move_to([(k-3.5)*1.55,.9,0]) for k in range(8)])
        self.play(FadeIn(bars),Create(Line([0,.3,0],[0,1.65,0],color=ORANGE)))
        group=VGroup(self.eq(r'd_k=(k-\tfrac{N_c-1}{2})P,\qquad h_{kj}=d_k+u_j-\bar u'),self.eq(r'T_{\rm coh}=N_cN_s/f_s=144\,\mu\mathrm{s}'),self.prose('8 chirps; 225 samples per chirp; 12.5 MHz sampling',27),self.prose('Reconstructed chirp period: 25.37 microseconds',25,MUTED)).arrange(DOWN,buff=.45).move_to(DOWN*1.25)
        self.play(FadeIn(group));self.wait(.7)

        self.start('Relation to the straight-line geometry')
        self.content([r'R(t)=\sqrt{(Vt-x_0)^2+y_0^2},\qquad X_\star=Vt_\star-x_0',r'r_0=\sqrt{X_\star^2+y_0^2},\qquad v_0=\frac{VX_\star}{r_0}',r'a_0=\frac{V^2y_0^2}{r_0^3}\geq0'],['V is along-track velocity; v0 is signed radial velocity.','Positive range curvature can coexist with decreasing range.'])

        self.start('The quadratic approximation is local')
        self.content([r'R(t)=r_0+v_0h+\tfrac12a_0h^2+\tfrac16j_0h^3+\cdots',r'\delta\phi\simeq-\frac{4\pi f_c}{c}\,\delta R'],['Check phase error over the entire acquired train.','A denser grid cannot repair Taylor truncation error.'])

        self.start('The implemented FMCW phase')
        self.content([r'q_{kj}=e^{i\phi_{kj}},\qquad \tau_{kj}=2R(t_{kj})/c',r'\phi_{kj}=2\pi\left[-(f_0+\gamma u_j)\tau_{kj}+\tfrac12\gamma\tau_{kj}^2\right]'],['The chirp slope mixes range and Doppler in fast-time frequency.','Both stored I/Q orientations are searched.'])

        self.start('Two FFT frequencies identify range and velocity')
        self.content([r'K_0=-\frac{2(f_0+\gamma\bar u)}{c}+\frac{4\gamma r_0}{c^2}',r'f_{\rm slow}=K_0v_0,\qquad f_{\rm fast}=K_0v_0-\frac{2\gamma r_0}{c}',r'r_0=\frac{c(f_{\rm slow}-f_{\rm fast})}{2\gamma},\qquad v_0=\frac{f_{\rm slow}}{K_0}'],['Physical search bounds resolve FFT aliases.'])

        self.start('Dominant sinusoids plus a phase correction')
        self.content([r'q=e^{i\phi_0}e^{2\pi i[f_{\rm fast}(u-\bar u)+f_{\rm slow}d]}e^{i\psi}',r'\psi=\phi-\phi_0-2\pi[f_{\rm fast}(u-\bar u)+f_{\rm slow}d]',r'z_{\rm corrected}=z\,e^{-i\psi}'],['The correction depends on r0, v0 and a0.','FFT searches the remaining dominant complex sinusoids.'])

        self.start('Curvature compensation restores coherent addition')
        circles=VGroup(*[Circle(radius=1.25,color=MUTED).move_to([x,.3,0]) for x in [-3,3]])
        before=VGroup(*[Arrow([-3,.3,0],[-3+1.15*np.cos(p),.3+1.15*np.sin(p),0],buff=0,color=ORANGE,stroke_width=2) for p in psi])
        after=Arrow([3,.3,0],[4.15,.3,0],buff=0,color=BLUE,stroke_width=7)
        labels=VGroup(self.prose('Before correction',28).move_to([-3,-1.4,0]),self.prose('After correction: aligned',28).move_to([3,-1.4,0]),self.prose('Synthetic phasors from eight chirps',25,MUTED).move_to([0,-2.6,0]))
        self.play(Create(circles),FadeIn(before),FadeIn(after),FadeIn(labels));self.wait(.7)

        self.start('The multi-chirp FFT search')
        self.content([],['Loop over acceleration planes and velocity correction groups.','Pre-apply the conjugate correction to each RX.','FFT across fast time, then across chirps.','Read bins for every explicit range and velocity node.','Add the four RX matched powers incoherently.','Search both I/Q orientations, then refine the peak.'])

        self.start('Human bounds in; grid spacing out')
        self.content([r'\epsilon=\arccos\sqrt{1-\ell},\qquad \max_{k,j}|\delta\psi_{kj}|\leq\epsilon',r'\Delta a_0\leq\frac{2\epsilon}{\max|\partial\psi/\partial a_0|}'],['Enter range, velocity and acceleration bounds and allowed loss.','FFT bin spacing sets range and velocity resolution.','Global phase bounds set acceleration spacing and correction groups.','Keep complete bounds; reject oversized grids rather than truncate.'])

        self.start('Matched powers and maximum projections')
        self.content([r'\rho=\frac{\sum_b|\langle q,z_b\rangle|^2}{\sum_bP_{n,b}\sum|q|^2}',r'M_{rv}(r_0,v_0)=\max_{a_0}\rho,\qquad M_{va}(v_0,a_0)=\max_{r_0}\rho'],['Trajectory search uses all four RX powers.','Beamforming refines relative phases after the trajectory fit.','A search peak alone does not establish a detection.'])

        self.start('Coherent time and noise bandwidth')
        self.content([r'T_{\rm coh}=\frac{N_cN_s}{f_s},\qquad B_{\rm analysis}=\frac1{T_{\rm coh}}',r'P_{n,b}^{\rm analysis}=P_{n,b}\frac{B_{\rm analysis}}{f_s}',r'N_c=8:\quad T_{\rm coh}=144\,\mu\mathrm{s},\quad B_{\rm analysis}=6.94\,\mathrm{kHz}'],["Estimate each RX's noise from quiet complex-mean-removed voltage.",'Use every intact sample in the selected background interval.'])

        self.start('The GUI: select background, select analysis, press Play','Schematic of actual controls, not a measurement screenshot. White display; full reset restores the view. CPU fallback is available when WebGPU is unavailable.')
        panel=Rectangle(width=12,height=2.5,color=MUTED).shift(UP*.4)
        bg=Rectangle(width=2,height=2.5,color=BLUE,fill_color=BLUE,fill_opacity=.14).move_to([-3.7,.4,0])
        scan=Rectangle(width=5.4,height=2.5,color=PURPLE,fill_color=PURPLE,fill_opacity=.10).move_to([1.4,.4,0])
        coh=Rectangle(width=.45,height=2.5,color=ORANGE,fill_color=ORANGE,fill_opacity=.2).move_to([-.7,.4,0])
        labels=VGroup(self.prose('Background',24,BLUE).move_to([-3.7,2,0]),self.prose('Analysis interval',24,PURPLE).move_to([1.4,2,0]),self.prose('Play',30,BLUE).move_to([-4,-1.5,0]),self.prose('Reset full view',27).move_to([0,-1.5,0]),self.prose('Download HDF5',27).move_to([4,-1.5,0]))
        notes=VGroup(self.prose('Move or resize blue and purple windows; drag empty space to zoom.',25),self.prose('Moving the yellow train analyzes that coherent interval.',25),self.prose('Results appear during the scan; HDF5 unlocks when it finishes.',25)).arrange(DOWN,buff=.27).move_to(DOWN*2.7)
        self.play(FadeIn(panel),FadeIn(bg),FadeIn(scan),FadeIn(coh),FadeIn(labels),FadeIn(notes));self.wait(.7)

        self.start('Reading the results and sharing the view')
        self.content([],['Maximum-projection maps: velocity versus range and acceleration.','Waveforms and fitted continuous range are line plots.','Scan histories are unconnected points: range, velocity, acceleration,','RX phases, RX and beamformed SNR, RCS and sphere diameters.','The URL preserves selections, bounds, timing and plot view.','Completed results download as HDF5.'],)
        # Keep the completed final slide visible; no outgoing fade.
