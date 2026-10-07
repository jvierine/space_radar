"""FYS-3000, 7 October 2026: radar power, thermal noise, coherent bandwidth.

Render: conda run -n base manim-slides render --disable_caching
        -r 1920,1080 --fps 30 radar_equation_noise.py RadarEquationNoise
The source footer is shown by default; set SHOW_PROVENANCE=0 to hide it.
"""
from __future__ import annotations

import os
import numpy as np
from manim import *
from manim_slides import Slide
from planck_to_ktb import Text

BG, FG, MUTED = '#07111F', '#F3F7FA', '#9DB0C3'
GOLD, SIGNAL, RECEIVE = '#F3D35A', '#74A9FF', '#45C2B1'
POWER, RED = '#F5A65B', '#FF6B6B'
SHOW_PROVENANCE = os.getenv('SHOW_PROVENANCE', '1') == '1'

RADAR_SOURCE = 'https://www.ll.mit.edu/sites/default/files/outreach/doc/2018-07/lecture%202.pdf'
NOISE_SOURCE = 'https://ocw.mit.edu/courses/6-013-electromagnetics-and-applications-spring-2009/d3be4ea78b036a6362230fb41780cf54_MIT6_013S09_notes.pdf'
FILTER_SOURCE = 'https://www.mathworks.com/help/signal/ref/enbw.html'

NOTES = {
    'Radar equation': 'Received available echo power, before receiver amplification. Monostatic, far field, point target, matched polarization. Gains are linear and evaluated toward the target. L>=1 contains losses not already included in antenna gain; avoid counting losses twice. For pulses Pt is on-pulse power; duty cycle and live signal time matter in later integration. Source: '+RADAR_SOURCE,
    'Transmit power and gain': 'Pt is power delivered to the transmit antenna. Gain includes antenna radiation efficiency. Gt multiplies directional flux relative to an isotropic antenna with the same delivered power. The first inverse-square factor is propagation from radar to target. Source: '+RADAR_SOURCE,
    'Radar cross section': 'Sigma=4*pi*R^2*S_back/S_inc in the far field, for the chosen direction and polarization. RCS is an equivalent scattering area, not necessarily the projected geometric area. It depends on wavelength, aspect, material and polarization. Source: '+RADAR_SOURCE,
    'Receiving aperture and losses': 'Ae is effective collection area, not necessarily physical antenna area. Ae=Gr*lambda^2/(4*pi) for the reciprocal matched receiving antenna. Wavelength and gain must refer to the same operating frequency. L is a linear power-loss factor. Two spreading factors produce R^-4. Source: '+RADAR_SOURCE,
    'Thermal noise power': 'Use available input-referred noise power and the Rayleigh-Jeans regime. Tsys includes antenna/environment and receiver-added noise, referred to the same reference plane as Pr. White thermal noise is assumed over the filter passband. This equation does not model clutter or RFI. Sources: '+RADAR_SOURCE+' ; '+NOISE_SOURCE,
    'Coherent integration': 'Complex I/Q convention: de-rotate by the predicted signal phase, then average amplitudes, not powers. A constant-amplitude signal is preserved by the normalized average. Unknown constant phase is allowed. Frequency offset, acceleration or range migration need matching or compensation. Integration must use samples containing the echo, with phase continuity; gaps require actual sample times. Source: '+NOISE_SOURCE,
    'Integration filter': 'After phase matching, h(t)=1/Tint for 0<=t<=Tint, zero otherwise, with unit DC gain. H(f)=exp(-i*pi*f*Tint)*sinc(f*Tint), where sinc(x)=sin(pi*x)/(pi*x). The plotted power response is normalized and uses two-sided complex-baseband frequency. First nulls are at +/-1/Tint. ENBW is the total area under |H|^2, not the null-to-null width or the 3-dB width. Source: '+FILTER_SOURCE,
    'Equivalent noise bandwidth': 'For complex I/Q, use the full two-sided baseband integral. Parseval gives integral |H|^2 df = integral |h|^2 dt = 1/Tint. In discrete time, B=fs*sum(|w|^2)/|sum(w)|^2; equal weights give fs/N=1/Tint. A Hann window gives approximately 1.5/Tint instead. A real one-sided low-pass convention may quote 1/(2*Tint); use a consistent PSD/bandwidth convention. Source: '+FILTER_SOURCE,
    'Ten seconds of coherent measurements': 'Example assumes Tsys=300 K, unit-gain rectangular coherent average, continuous signal and wide input passband. kB is exact in SI. Noise=4.141947e-22 W at 0.1 Hz. Changing integration from 1 s to 10 s reduces output noise power by ten, while leaving a matched constant signal unchanged. Source: '+NOISE_SOURCE+' ; '+FILTER_SOURCE,
    'Radar signal to noise ratio': 'The displayed SNR assumes a constant received power throughout the integration and perfect phase matching. For varying received power the matched-filter SNR uses received signal energy divided by kB*Tsys, with appropriate weighting. Longer recording alone does not improve the SNR of an echo that has already ended. Sources: '+RADAR_SOURCE+' ; '+FILTER_SOURCE,
}


class RadarEquationNoise(Slide):
    def setup(self):
        self._slide_started = False
        self._clear_pending = False
        self.camera.background_color = ManimColor(BG)
        self.provenance = Text('Source: radar_equation_noise.py', font_size=14, color=MUTED).to_corner(DR, buff=.15)
        if SHOW_PROVENANCE:
            self.add(self.provenance)

    def start_slide(self, name):
        if self._slide_started:
            self.next_slide(name, notes=NOTES[name])
            self._flush_pending_clear()
        else:
            self.next_slide(name, notes=NOTES[name])
            self._slide_started = True

    def clear_slide(self):
        self._clear_pending = True

    def _flush_pending_clear(self):
        if self._clear_pending:
            self._clear_pending = False
            keep = {self.provenance} if SHOW_PROVENANCE else set()
            removable = [mob for mob in self.mobjects if mob not in keep]
            if removable:
                self.play(*[FadeOut(mob) for mob in removable], run_time=.45)

    def title(self, text):
        result = Text(text, font_size=43, weight=SEMIBOLD, color=FG)
        if result.width > 13:
            result.scale_to_fit_width(13)
        return result.to_edge(UP, buff=.32)

    def prose(self, text, size=29, color=FG):
        result = Text(text, font_size=size, color=color)
        if result.width > 12.6:
            result.scale_to_fit_width(12.6)
        return result

    def equation(self, *tex, size=50, color=FG):
        result = MathTex(*tex, font_size=size, color=color)
        if result.width > 12.8:
            result.scale_to_fit_width(12.8)
        return result

    def finish(self):
        self.wait(.8)
        self.clear_slide()

    def geometry(self):
        radar = VGroup(Line([-5,-.8,0],[-5,.2,0],color=SIGNAL,stroke_width=5),
                       Arc(radius=.45,start_angle=-PI/2,angle=PI,color=SIGNAL).move_to([-4.75,.15,0]))
        target = Circle(radius=.35,color=GOLD,fill_color=GOLD,fill_opacity=.35).move_to([4.5,.1,0])
        outgoing = Arrow([-4.1,.4,0],[3.9,.4,0],buff=0,color=POWER)
        incoming = Arrow([3.9,-.25,0],[-4.1,-.25,0],buff=0,color=RECEIVE)
        labels = VGroup(self.prose('Radar',25,SIGNAL).move_to([-5,-1.15,0]),
                        self.prose('Target',25,GOLD).move_to([4.5,-1.15,0]),
                        self.equation('R',size=32,color=MUTED).move_to([0,.85,0]))
        return VGroup(radar,target,outgoing,incoming,labels)

    def construct(self):
        self.radar_equation()
        self.transmit()
        self.cross_section()
        self.aperture()
        self.noise()
        self.coherent()
        self.filter_response()
        self.bandwidth()
        self.ten_seconds()
        self.snr()
        # Leave the final completed slide visible. No outgoing cleanup.

    def radar_equation(self):
        self.start_slide('Radar equation')
        self.play(FadeIn(self.title('The radar equation')))
        eq = self.equation(r'P_r=\frac{P_tG_tG_r\lambda^2\sigma}{(4\pi)^3R^4L}',size=64).move_to(UP*1.55)
        self.play(Write(eq),run_time=1.6)
        diagram = self.geometry().shift(DOWN*.55)
        self.play(LaggedStart(*[FadeIn(part) for part in diagram],lag_ratio=.15),run_time=1.5)
        text = self.prose('Received echo power = illumination × scattering × collection / losses',29).move_to(DOWN*2.45)
        self.play(FadeIn(text))
        note = self.prose('Monostatic radar: transmitter and receiver at the same location',23,MUTED).move_to(DOWN*3.15)
        self.play(FadeIn(note))
        self.finish()

    def transmit(self):
        self.start_slide('Transmit power and gain')
        self.play(FadeIn(self.title('Illumination of the target')))
        eq = self.equation(r'S_{\rm inc}=\frac{P_tG_t}{4\pi R^2}',size=62,color=POWER).move_to(UP*1.6)
        self.play(Write(eq))
        rows = VGroup(
            self.prose('Transmit power: watts delivered to the antenna',31),
            self.prose('Transmit gain: concentration toward the target',31),
            self.prose('Range in metres: power spreads over an area proportional to range squared',29),
        ).arrange(DOWN,aligned_edge=LEFT,buff=.52).move_to(DOWN*.3)
        symbols = VGroup(*[self.equation(s,size=40,color=POWER) for s in [r'P_t',r'G_t',r'R']])
        for symbol,row in zip(symbols,rows):
            symbol.next_to(row,LEFT,buff=.45)
        group = VGroup(symbols,rows)
        if group.width > 12.6: group.scale_to_fit_width(12.6)
        group.move_to(DOWN*.35)
        for symbol,row in zip(symbols,rows): self.play(FadeIn(symbol),FadeIn(row),run_time=.65)
        self.play(FadeIn(self.equation(r'[S_{\rm inc}]=\mathrm{W\,m^{-2}},\qquad G_t\text{ is a linear power ratio}',size=29,color=MUTED).move_to(DOWN*2.7)))
        self.finish()

    def cross_section(self):
        self.start_slide('Radar cross section')
        self.play(FadeIn(self.title('Radar cross section')))
        definition = self.equation(r'\sigma=4\pi R^2\frac{S_{\rm back}}{S_{\rm inc}}',size=56,color=GOLD).move_to(UP*1.65)
        self.play(Write(definition))
        self.play(FadeIn(self.prose('An equivalent area describing scattering back toward the radar',31).move_to(UP*.3)))
        self.play(FadeIn(self.prose('Units: square metres. It can differ greatly from geometric area.',29).move_to(DOWN*.4)))
        self.play(FadeIn(self.prose('Depends on wavelength, shape, orientation, material and polarization',27,MUTED).move_to(DOWN*1.15)))
        return_eq = self.equation(r'S_{\rm back}=S_{\rm inc}\frac{\sigma}{4\pi R^2}',size=49,color=RECEIVE).move_to(DOWN*2.35)
        self.play(Write(return_eq))
        self.finish()

    def aperture(self):
        self.start_slide('Receiving aperture and losses')
        self.play(FadeIn(self.title('Collection by the receiving antenna')))
        eq = self.equation(r'A_e=\frac{G_r\lambda^2}{4\pi}',size=53,color=RECEIVE).move_to(UP*1.85)
        self.play(Write(eq))
        self.play(FadeIn(self.prose('Effective aperture: the area that collects the returning flux',29).move_to(UP*.75)))
        definitions = VGroup(self.prose('Receive gain',27),self.equation(r'G_r',size=33,color=RECEIVE),
                             self.prose('and wavelength',27),self.equation(r'\lambda=c/f\quad[\mathrm m]',size=33,color=RECEIVE)).arrange(RIGHT,buff=.22).move_to(DOWN*.05)
        self.play(FadeIn(definitions))
        product = self.equation(r'P_r=\underbrace{\frac{P_tG_t}{4\pi R^2}}_{\text{illumination}}\quad'
                                r'\underbrace{\frac{\sigma}{4\pi R^2}}_{\text{scattering}}\quad'
                                r'\underbrace{\frac{G_r\lambda^2}{4\pi}}_{\text{collection}}\quad\frac{1}{L}',size=43).move_to(DOWN*1.25)
        self.play(Write(product),run_time=1.8)
        note = self.prose('L is a linear loss factor: L = 1 is ideal; L = 2 halves the echo power',25,MUTED).move_to(DOWN*2.4)
        self.play(FadeIn(note))
        scaling = self.equation(r'P_r\propto R^{-4}:\quad 2R\ \Longrightarrow\ P_r/16',size=37,color=GOLD).move_to(DOWN*3)
        self.play(Write(scaling))
        self.finish()

    def noise(self):
        self.start_slide('Thermal noise power')
        self.play(FadeIn(self.title('Thermal noise power')))
        self.play(Write(self.equation(r'P_n=k_B T_{\rm sys}B_{\rm eff}',size=68,color=POWER).move_to(UP*1.65)))
        rows = [
            (r'k_B=1.380649\times10^{-23}\ \mathrm{J\,K^{-1}}','Boltzmann constant'),
            (r'T_{\rm sys}\ [\mathrm K]','Antenna noise and receiver-added noise'),
            (r'B_{\rm eff}\ [\mathrm{Hz}]','Equivalent noise bandwidth of the measurement filter'),
        ]
        for i,(symbol,description) in enumerate(rows):
            y=.25-i*.85
            label=self.equation(symbol,size=31,color=POWER).move_to([-3.2,y,0])
            text=self.prose(description,27).move_to([2.65,y,0])
            if text.width>6: text.scale_to_fit_width(6)
            self.play(FadeIn(label),FadeIn(text),run_time=.7)
        self.play(FadeIn(self.equation(r'N_0=k_BT_{\rm sys}\quad[\mathrm{W\,Hz^{-1}}]',size=34,color=MUTED).move_to(DOWN*2.7)))
        self.finish()

    def coherent(self):
        self.start_slide('Coherent integration')
        self.play(FadeIn(self.title('Coherent integration')))
        model=self.equation(r'z(t)=A e^{i\phi(t)}+n(t)',size=48).move_to(UP*1.8)
        self.play(Write(model))
        estimate=self.equation(r'\widehat A=\frac{1}{T_{\rm int}}\int_0^{T_{\rm int}}z(t)e^{-i\phi(t)}\,dt',size=52,color=SIGNAL).move_to(UP*.4)
        self.play(Write(estimate))
        self.play(FadeIn(self.prose('Match the expected phase, then average complex voltages',32).move_to(DOWN*.7)))
        self.play(FadeIn(self.prose('The aligned signal adds coherently; random noise averages down',29).move_to(DOWN*1.45)))
        self.play(FadeIn(self.prose('Motion and Doppler must be accounted for during the integration',26,MUTED).move_to(DOWN*2.45)))
        self.finish()

    def filter_response(self):
        self.start_slide('Integration filter')
        self.play(FadeIn(self.title('Longer integration gives a narrower filter')))
        response=self.equation(r'|H(f)|^2=\left[\frac{\sin(\pi fT_{\rm int})}{\pi fT_{\rm int}}\right]^2',size=38).move_to(UP*2.25)
        self.play(Write(response))
        ax=Axes(x_range=[-1.1,1.1,.5],y_range=[0,1.05,.5],x_length=10.2,y_length=2.7,
                tips=False,axis_config={'color':MUTED},x_axis_config={'include_numbers':True},y_axis_config={'include_numbers':True}).move_to(DOWN*.15)
        xlab=self.prose('Frequency offset after phase matching (Hz)',24,MUTED).next_to(ax,DOWN,buff=.28)
        ylab=self.equation(r'|H|^2',size=26,color=MUTED).next_to(ax,LEFT,buff=.25)
        curve1=ax.plot(lambda f:np.sinc(f)**2,x_range=[-1.1,1.1,.005],color=SIGNAL,stroke_width=4,use_smoothing=False)
        curve10=ax.plot(lambda f:np.sinc(10*f)**2,x_range=[-1.1,1.1,.002],color=GOLD,stroke_width=4,use_smoothing=False)
        label1=self.equation(r'T_{\rm int}=1\ \mathrm s',size=28,color=SIGNAL).move_to([-3.5,1.15,0])
        label10=self.equation(r'T_{\rm int}=10\ \mathrm s',size=28,color=GOLD).move_to([3.5,1.15,0])
        self.play(Create(ax),FadeIn(xlab),FadeIn(ylab))
        self.play(Create(curve1),FadeIn(label1),run_time=1.2)
        self.play(Create(curve10),FadeIn(label10),run_time=1.2)
        self.play(FadeIn(self.prose('A ten-second average rejects frequency offsets much more strongly',27).move_to(DOWN*2.7)))
        self.finish()

    def bandwidth(self):
        self.start_slide('Equivalent noise bandwidth')
        self.play(FadeIn(self.title('Where the noise bandwidth comes from')))
        definition=self.equation(r'B_{\rm eff}=\frac{\int_{-\infty}^{\infty}|H(f)|^2\,df}{|H(0)|^2}',size=51).move_to(UP*1.5)
        self.play(Write(definition))
        self.play(FadeIn(self.prose('The area under the normalized filter power response',30).move_to(UP*.25)))
        rect=self.equation(r'h(t)=\frac{1}{T_{\rm int}}\quad(0\le t\le T_{\rm int})',size=41,color=SIGNAL).move_to(DOWN*.75)
        result=self.equation(r'B_{\rm eff}=\int |h(t)|^2\,dt=\frac{1}{T_{\rm int}}',size=51,color=GOLD).move_to(DOWN*1.8)
        self.play(Write(rect)); self.play(Write(result))
        self.play(FadeIn(self.prose('Rectangular coherent average, complex I/Q convention',24,MUTED).move_to(DOWN*2.95)))
        self.finish()

    def ten_seconds(self):
        self.start_slide('Ten seconds of coherent measurements')
        self.play(FadeIn(self.title('Ten seconds of coherent measurements')))
        self.play(Write(self.equation(r'T_{\rm int}=10\ \mathrm s\quad\Longrightarrow\quad B_{\rm eff}=0.1\ \mathrm{Hz}',size=57,color=GOLD).move_to(UP*1.75)))
        self.play(FadeIn(self.prose('Example: system noise temperature of 300 K',32).move_to(UP*.4)))
        numeric=self.equation(r'P_n=(1.380649\times10^{-23})(300)(0.1)',r'=4.14\times10^{-22}\ \mathrm W',size=45,color=POWER).arrange(DOWN,buff=.28).move_to(DOWN*.85)
        self.play(Write(numeric))
        self.play(FadeIn(self.prose('Ten times longer integration gives ten times less noise power',29).move_to(DOWN*2.3)))
        self.play(FadeIn(self.prose('The signal must be present and correctly phase-matched throughout',25,MUTED).move_to(DOWN*3)))
        self.finish()

    def snr(self):
        self.start_slide('Radar signal to noise ratio')
        self.play(FadeIn(self.title('Radar signal-to-noise ratio')))
        equation=self.equation(r'\mathrm{SNR}=\frac{P_r}{k_BT_{\rm sys}B_{\rm eff}}',size=60,color=SIGNAL).move_to(UP*1.8)
        self.play(Write(equation))
        combined=self.equation(r'\mathrm{SNR}=\frac{P_tG_tG_r\lambda^2\sigma\,T_{\rm int}}{(4\pi)^3R^4L\,k_BT_{\rm sys}}',size=55).move_to(UP*.05)
        self.play(Write(combined),run_time=1.5)
        self.play(FadeIn(self.prose('Echo strength is set by the radar, the target and the range',31).move_to(DOWN*1.5)))
        self.play(FadeIn(self.prose('A longer coherent measurement reduces the noise bandwidth',30).move_to(DOWN*2.25)))
        self.play(FadeIn(self.prose('Constant echo power, rectangular weighting, correct phase matching',23,MUTED).move_to(DOWN*3.05)))
        self.wait(1)
