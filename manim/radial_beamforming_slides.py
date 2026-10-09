"""Beamforming slides by Juha Vierinen.
PCB layout: TI SPRUIM4B figure 10, https://www.ti.com/lit/ug/spruim4b/spruim4b.pdf.
Geometry schematic; plane-wave phases are synthetic, not measurement results.
"""
import numpy as np
from manim import *
BLUE, ORANGE, PURPLE, MUTED = '#176BB0','#B55A00','#7939A8','#556575'
COLORS=[BLUE,ORANGE,'#008777',PURPLE]


def beamforming_slides(s):
    s.start('The four receivers are four PCB antenna columns',
        'TI SPRUIM4B figure 10. Four RX columns with three series-fed patches each, flanked by dummy columns. Adjacent RX pitch lambda/2. Three TX columns at right, horizontal pitch lambda; middle TX raised lambda/2. Simplified PCB layout, not fabrication dimensions. TI labels RX1-RX4; GUI uses RX0-RX3. Capture-to-physical channel order needs verification before angle calibration.')
    board=RoundedRectangle(width=11.4,height=3.3,corner_radius=.15,stroke_color='#993B41',fill_color='#F8E9EA',fill_opacity=1).move_to(UP*.2)
    s.play(FadeIn(board))
    patches=VGroup();labels=VGroup()
    columns=[(-4.6,'Dummy',MUTED,0),(-3.6,'RX0',BLUE,0),(-2.6,'RX1',ORANGE,0),(-1.6,'RX2',COLORS[2],0),(-.6,'RX3',PURPLE,0),(.4,'Dummy',MUTED,0),(2.4,'TX1',MUTED,0),(4.4,'TX2',MUTED,.5),(6.4,'TX3',MUTED,0)]
    for x,label,color,dy in columns:
        x=(x-.9)*.82
        for y in [-.55,0,.55]:
            patches.add(Rectangle(width=.37,height=.22,stroke_color=color,fill_color=color,fill_opacity=.8).move_to([x,y+dy+.3,0]))
        patches.add(Line([x,-1,0],[x,.9+dy,0],color=color,stroke_width=2))
        labels.add(s.prose(label,20,color).move_to([x,1.45+dy*.25,0]))
    s.play(FadeIn(patches),FadeIn(labels))
    s.play(Write(s.eq(r'd=\lambda/2\approx1.95\ \mathrm{mm}\quad(77\ \mathrm{GHz})',35).move_to(DOWN*1.85)))
    text=VGroup(s.prose('Each RX column feeds one synchronized complex voltage channel.',28),s.prose('Three patches per column; dummy columns are not recorded receivers.',25,MUTED),s.prose('xWR1843BOOST PCB schematic, based on TI figure 10.',21,MUTED)).arrange(DOWN,buff=.24).move_to(DOWN*2.85)
    s.play(FadeIn(text));s.wait(1)

    s.start('An incoming wave gives each antenna a different phase',
        'Top view of receive phase centers, far-field plane wave theta=30 degrees from broadside. Convention exp(+i omega t), incoming propagation toward negative x and y: phase +2pi x sin(theta)/lambda. Adjacent phase pi sin(theta), here 90 degrees. Conjugated IQ reverses sign. This is receive-only path difference, no monostatic factor of two. Electronic phase offsets also contribute.')
    xs=[-3.,-1.,1.,3.];base=-.6;theta=PI/6
    receivers=VGroup(*[Dot([x,base,0],radius=.13,color=COLORS[b]) for b,x in enumerate(xs)])
    labels=VGroup(*[s.prose(f'RX{b}',25,COLORS[b]).move_to([x,base-.4,0]) for b,x in enumerate(xs)])
    def wavefront(y):
        slope=np.tan(theta)
        x0=max(-3.8,(y-2.3)/slope);x1=min(3.8,(y+.6)/slope)
        return Line([x0,y-slope*x0,0],[x1,y-slope*x1,0],color=BLUE,stroke_opacity=.55)
    offsets=[.1,.85,1.6]
    fronts=VGroup(*[wavefront(y) for y in offsets])
    broadside=DashedLine([0,base,0],[0,2.55,0],color=MUTED)
    incoming=Arrow([1.6,2.2,0],[.6,.47,0],color=ORANGE,buff=0)
    angle=Arc(radius=.8,start_angle=PI/2-theta,angle=theta,arc_center=[0,base,0],color=ORANGE)
    s.play(FadeIn(receivers),FadeIn(labels),FadeIn(broadside),FadeIn(fronts),GrowArrow(incoming),Create(angle))
    s.play(FadeIn(s.eq(r'\theta',26).move_to([.35,.45,0])),FadeIn(s.prose('Incoming wave',23,ORANGE).move_to([3,2.3,0])))
    s.play(*[Transform(front,wavefront(y-.6/np.cos(theta))) for front,y in zip(fronts,offsets)],run_time=1.8,rate_func=linear)
    s.play(Write(s.eq(r'\Delta\ell=d\sin\theta,\qquad\Delta\phi=\frac{2\pi d}{\lambda}\sin\theta',37).move_to(DOWN*1.65)))
    definitions=VGroup(s.prose('theta: arrival angle from the array normal (broadside).',25),s.prose('d: spacing. Delta ell: extra receive path. Delta phi: phase difference.',23),s.prose('Half-wavelength spacing: 30 degrees arrival → 90 degrees between RXs.',25,BLUE)).arrange(DOWN,buff=.25).move_to(DOWN*2.85)
    s.play(FadeIn(definitions));s.wait(1)

    s.start('Rotate the matched outputs, then add complex amplitudes',
        'At one common best trajectory sb=sum q* zb. Preserve complex numbers; powers alone cannot beamform. Synthetic equal outputs phases 0,90,180,270 degrees; rotate alpha=-phase relative to RX0. These circles show four matched-filter outputs, not chirps. B=(1/2) sum exp(i alpha_b) sb has unit-norm weights. Linearity allows beamforming before or after the same temporal matched filter. The displayed analytic alignment maximizes power for equal independent noise; actual GUI maximizes noise-normalized score.')
    centres=[np.array([x,.85,0]) for x in [-4.5,-1.5,1.5,4.5]]
    circles=VGroup(*[Circle(radius=.75,color=MUTED).move_to(c) for c in centres])
    arrows=VGroup(*[Arrow(c,c+.68*np.array([np.cos(b*PI/2),np.sin(b*PI/2),0]),buff=0,color=COLORS[b]) for b,c in enumerate(centres)])
    labels=VGroup(*[s.prose(f'RX{b}',26,COLORS[b]).move_to(c+UP*1.15) for b,c in enumerate(centres)])
    s.play(FadeIn(circles),FadeIn(labels),*[GrowArrow(a) for a in arrows])
    s.play(Write(s.eq(r's_b=\sum_jq_j^*z_{b,j}\quad\text{at the same best }(r_0,v_0,a_0)',35).move_to(UP*2.65)))
    caption=s.prose('Before rotation: the four complex outputs cancel in this example.',27,ORANGE).move_to(DOWN*.5)
    s.play(FadeIn(caption));s.wait(.5)
    s.play(*[Rotate(a,angle=-b*PI/2,about_point=c) for b,(a,c) in enumerate(zip(arrows,centres))],run_time=2)
    s.play(Transform(caption,s.prose('After rotation: all four echo amplitudes add constructively.',27,BLUE).move_to(DOWN*.5)))
    s.play(Write(s.eq(r'B=\frac12\sum_{b=0}^3e^{i\alpha_b}s_b,\qquad\alpha_b=-\arg(s_b)+\arg(s_0)',38).move_to(DOWN*1.55)))
    terms=VGroup(s.prose('s(b): complex matched output. alpha(b): applied phase rotation.',24),s.prose('B: beamformed output. One-half gives unit-norm four-channel weights.',24),s.prose('Equal aligned channels with independent noise give four times the single-RX SNR.',24,BLUE)).arrange(DOWN,buff=.23).move_to(DOWN*2.9)
    s.play(FadeIn(terms));s.wait(1)

    s.start('Search the receiver phases and propagate the background noise',
        'Implementation fixes RX0=0, coarse three-phase grid then Nelder-Mead refinement of observed beam power over predicted noise. C_ab=E[n_a conjugate(n_b)] is matched-output noise covariance, quiet sample covariance times template energy. c_b=exp(i alpha_b)/2; beam noise=sum c_a conjugate(c_b) C_ab. Equal independent noise: noise unchanged, signal power fourfold. Unequal channels and correlated noise change gain. General covariance-optimal amplitude weights are not implemented. Fitted phases are not calibrated arrival angles.')
    s.content([r'P_{n,\mathrm{beam}}=\sum_{a,b}c_a c_b^*C_{ab},\qquad c_b=\tfrac12e^{i\alpha_b}',r'\rho_{\mathrm{beam}}=\frac{|B|^2}{P_{n,\mathrm{beam}}}',r'\text{equal independent RXs:}\quad G_{\mathrm{SNR}}=4\approx6.02\ \mathrm{dB}'],[
        'C: matched-output noise covariance, estimated from quiet background.',
        'Fix RX0 at zero; search three relative phases, then refine the best combination.',
        'Carry each channel noise and cross-channel correlations through the same weights.',
        'Receiver phases include geometry and electronic offsets.',
        'Unequal channels can benefit from amplitude weighting; this GUI fits phases only.',
    ])
