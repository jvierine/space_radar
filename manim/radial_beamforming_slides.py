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
    definitions=VGroup(s.texrow(r'$\theta$: arrival angle from the array normal (broadside).',25),s.texrow(r'$d$: spacing; $\Delta\ell$: extra receive path; $\Delta\phi$: phase difference.',23),s.texrow(r'With $d=\lambda/2$: $\theta=30^\circ$ gives $\Delta\phi=90^\circ$ between RXs.',25).set_color(BLUE)).arrange(DOWN,buff=.25).move_to(DOWN*2.85)
    s.play(FadeIn(definitions));s.wait(1)

    s.start('Rotate the matched outputs, then add complex amplitudes',
        'At one common best trajectory sb=sum q* zb. Preserve complex numbers; powers alone cannot beamform. Synthetic equal outputs phases 0,90,180,270 degrees; rotate alpha=-phase relative to RX0. These circles show four matched-filter outputs, not chirps. s_{\rm BF}=(1/2) sum exp(i alpha_b) sb has unit-norm weights. Linearity allows beamforming before or after the same temporal matched filter. The displayed analytic alignment maximizes power for equal independent noise; This equal-noise illustration is the special case of the joint noise-weighted fit shown next.')
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
    s.play(Write(s.eq(r's_{\rm BF}=\frac12\sum_{b=0}^3e^{i\alpha_b}s_b,\qquad\alpha_b=-\arg(s_b)+\arg(s_0)',38).move_to(DOWN*1.55)))
    terms=VGroup(s.texrow(r'$s_b$: complex matched output; $\alpha_b$: applied phase rotation.',24),s.texrow(r'$s_{\mathrm{BF}}$: beamformed output; $1/2$ gives unit-norm four-channel weights.',24),s.prose('Equal aligned channels with independent noise give four times the single-RX SNR.',24,BLUE)).arrange(DOWN,buff=.23).move_to(DOWN*2.9)
    s.play(FadeIn(terms));s.wait(1)

    s.start('One joint fit refines motion and the complex receiver amplitudes',
        'Joint generalized least squares. Sample vector z_j has four complex receivers, z_j=q_j(theta) A+n_j. C is per-sample receive covariance estimated in the quiet residual; temporally white acquired noise assumed. Minimize sum (z_j-q_j A)^H C^-1 (z_j-q_j A). For common q and fixed C, derivative yields A=s/Q, s=sum q* z and Q=sum|q|^2. The profiled objective equals constant minus s^H C^-1 s/Q. Bounded Nelder-Mead adjusts theta; amplitudes are solved exactly each trial. This includes motion refinement and receiver phase/amplitude fitting in the same objective, no subsequent phase-only simplex. Search peaks retain selection bias.')
    model=MathTex(r'\mathbf z_j=',r'q_j(r_0,v_0,a_0)',r'\mathbf A',r'+\mathbf n_j',font_size=44,color='#202830').move_to(UP*1.8)
    s.play(Write(model))
    s.callout(model[1],'Predicted echo, including motion',[-2.8,2.85,0],ORANGE)
    s.callout(model[2],'Four complex RX amplitudes',[3.7,2.85,0],BLUE)
    objective=s.eq(r'J(\theta,\mathbf A)=\sum_j(\mathbf z_j-q_j\mathbf A)^H C^{-1}(\mathbf z_j-q_j\mathbf A)',38).move_to(UP*.2)
    s.play(Write(objective))
    s.play(Write(s.eq(r'\mathbf s=\sum_jq_j^*\mathbf z_j,\quad Q=\sum_j|q_j|^2,\quad\widehat{\mathbf A}=\mathbf s/Q',37).move_to(DOWN*.9)))
    rows=VGroup(s.texrow(r'$C$: background receive-noise covariance; $H$: conjugate transpose.',25),
        s.prose('At each trial: solve amplitudes, weight residuals, evaluate least squares.',26),
        s.prose('Nelder–Mead adjusts range, radial velocity and radial acceleration together.',25,PURPLE),
        s.prose('Those same fitted amplitudes supply the beam phases and amplitude weights.',25,BLUE)).arrange(DOWN,buff=.25).move_to(DOWN*2.55)
    s.play(FadeIn(rows));s.wait(1)

    s.start('The same fit supplies the noise-weighted coherent beam',
        'Optimal linear combining for a known or fitted amplitude vector A under receive covariance C: w proportional C^-1 A, beam s_{\rm BF}=w^H s, matched noise Q w^H C w. Unit Euclidean normalization is arbitrary for SNR. Estimated weights use same event amplitudes; score s^H (Q C)^-1 s has noise-only expectation tr(Creg^-1 C), equal number of fitted complex RX amplitudes without regularization (four). Remove this amplitude-fitting noise bias when estimating excess signal SNR. Further trajectory search selection bias remains. RCS gain propagates actual weights and covariance assuming equal calibrated physical antenna responses; not universally four. Fitted phases are not calibrated arrival angles.')
    s.content([r'\mathbf w\propto C^{-1}\widehat{\mathbf A},\qquad s_{\rm BF}=\mathbf w^H\mathbf s',
        r'P_{n,\rm BF}=Q\mathbf w^H C\mathbf w,\qquad\rho_{\rm GLS}=\frac{\mathbf s^H C^{-1}\mathbf s}{Q}',
        r'\widehat S=\max(\rho_{\rm GLS}-\nu,0),\qquad\nu=\mathrm{tr}(C_{\rm reg}^{-1}C)\simeq4'],[
        r'$\mathbf w$: complex combining weights; $s_{\mathrm{BF}}$: beam voltage; $Q$: template energy.',
        'Inverse covariance downweights noisy channels and accounts for noise correlations.',
        r'$\nu$: noise contribution from fitting four complex receiver amplitudes.',
        'Independent equal-noise, equal-signal RXs give the familiar fourfold SNR gain.',
        'Search selection can still bias weak peaks; fitted phases include electronic offsets.',
    ])
