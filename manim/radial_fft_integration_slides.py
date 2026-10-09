"""Full-vector versus separable FFT integration, by Juha Vierinen.
The cost comparison is analytical; no equivalent zero-filled FFT benchmark is claimed.
"""
from manim import *
BLUE, ORANGE, PURPLE, MUTED = '#176BB0', '#B55A00', '#7939A8', '#556575'


def integration_slides(s):
    s.start('Two FFTs still integrate the entire complex chirp train',
        'For corrected data x_kj=z_kj c_kj conjugate, the separable template leaves fast and slow sinusoids. First F_k sums acquired fast-time samples with their complex phase, then S sums the complex F_k over chirps. Expanding yields the same double sum as a full-vector inner product, up to the common constant phase. FFT evaluation exact at its bin frequencies for that correction. Group sharing and off-bin lookup introduce mismatch in the search; full-template checks refine candidates. No absolute-square or incoherent integration between the FFT stages. User long-vector inner product is correct.')
    s.content([
        r'F_k(f_{\rm fast})=\sum_jx_{k,j}e^{-2\pi i f_{\rm fast}u_j}',
        r'S=\sum_kF_k(f_{\rm fast})e^{-2\pi i f_{\rm slow}d_k}',
        r'S=\sum_{k,j}x_{k,j}e^{-2\pi i(f_{\rm fast}u_j+f_{\rm slow}d_k)}',
    ],[
        'x: corrected complex samples. k: chirp. j: sample within that chirp.',
        'u: within-chirp time. d: chirp start offset. F: complex first-FFT output.',
        'Keep phase between stages; take squared magnitude only after the second FFT.',
        'Same matched sum and SNR as the corresponding long-vector inner product.',
        'A correct echo template maximizes SNR under white temporal noise.',
        'Colored temporal noise requires whitening; our current filter assumes white noise.',
    ])
    s.start('First FFT: every bin. Second FFT: the columns we need.',
        'Schematic, not measured data. Rows are chirps; columns are within-chirp frequency bins. Every row FFT computes every bin and retains complex values. The physical range/velocity grid maps to fast-frequency bins modulo the sampling rate. Deduplicate that list, then FFT every requested column across all acquired chirps. No selection by peak height. Aliasing is included in the modulo mapping; it is not why arbitrary columns can be discarded.')
    cells=VGroup()
    selected=[2,3,4,5]
    for k in range(4):
        for j in range(10):
            cell=Rectangle(width=.62,height=.5,stroke_color=MUTED,
                fill_color=BLUE if j in selected else '#DEE5EB',fill_opacity=.45).move_to([-4.7+j*.67,1.6-k*.6,0])
            cells.add(cell)
    row_labels=VGroup(*[s.prose(f'Chirp {k+1}',20).move_to([-6,1.6-k*.6,0]) for k in range(4)])
    s.play(FadeIn(cells),FadeIn(row_labels))
    arrows=VGroup(*[Arrow([-5.1,1.6-k*.6,0],[1.7,1.6-k*.6,0],color=ORANGE,buff=0,stroke_width=3) for k in range(4)])
    s.play(LaggedStart(*[GrowArrow(a) for a in arrows],lag_ratio=.25))
    s.play(FadeOut(arrows))
    s.play(FadeIn(s.prose('All frequency bins',25).move_to([-1.7,2.25,0])))
    vertical=VGroup(*[Arrow([-4.7+j*.67,1.95,0],[-4.7+j*.67,-.65,0],color=BLUE,buff=0) for j in selected])
    s.play(*[GrowArrow(a) for a in vertical])
    box=RoundedRectangle(width=3.3,height=1.4,stroke_color=BLUE).move_to([4,.7,0])
    label=VGroup(s.prose('FFT across chirps',25,BLUE),s.prose('Keep complex phase',22)).arrange(DOWN,buff=.25).move_to(box)
    s.play(FadeIn(box),FadeIn(label),GrowArrow(Arrow([1.9,.7,0],[2.3,.7,0],buff=0,color=BLUE)))
    rows=VGroup(s.prose('Blue columns: every bin requested by the range–velocity grid.',27,BLUE),
        s.prose('Grey columns: computed by the first FFT, unused by this search.',25),
        s.prose('Selection follows the physical search bounds, never signal strength.',25),
        s.prose('All chirps contribute coherently to every selected column.',26)).arrange(DOWN,buff=.28).move_to(DOWN*2.15)
    s.play(FadeIn(rows));s.wait(1)

    s.start('One beat frequency cannot separate range and velocity',
        'FMCW range-Doppler coupling. In our IQ convention, corrected fast frequency = slow frequency - 2 gamma r0/c. Approximate slow frequency = -2 fc v0/c, fc=f0+gamma mean fast time. Exact implemented slow expression has the small delay-squared range term: (-2fc/c+4gamma r0/c^2)*v0. Graphic is a synthetic constant-fast-frequency contour, approximately linear. Along the contour range compensates Doppler. These relations describe the corrected separable template, not a claim that an uncorrected accelerating signal is an exact tone.')
    axes=Axes(x_range=[.7,1.3,.2],y_range=[0,600,200],x_length=6.1,y_length=3.1,
        axis_config={'color':MUTED,'include_tip':False}).move_to([-2.7,.35,0])
    line=Line(axes.c2p(1.2,0),axes.c2p(.732,600),color=ORANGE,stroke_width=5)
    labels=VGroup(s.prose('Range r0 (m)',23).next_to(axes,DOWN),s.prose('Velocity v0 (m/s)',23).next_to(axes,UP))
    s.play(Create(axes),FadeIn(labels),Create(line))
    dot=Dot(axes.c2p(1.2,0),color=ORANGE)
    s.play(FadeIn(dot));s.play(MoveAlongPath(dot,line),run_time=2.5)
    equations=VGroup(s.eq(r'f_{\rm fast}\simeq-\frac{2\gamma r_0}{c}-\frac{2f_c v_0}{c}',32),
        s.prose('One measured frequency:',25,ORANGE),s.prose('many range–velocity pairs.',25,ORANGE)).arrange(DOWN,buff=.4).move_to([3.65,.5,0])
    s.play(FadeIn(equations))
    rows=VGroup(s.prose('Orange line: different targets with the same within-chirp beat frequency.',25),
        s.prose('gamma: chirp slope. fc: carrier near the sampled chirp midpoint.',24),
        s.prose('c: speed of light. The sign follows our complex I/Q convention.',24)).arrange(DOWN,buff=.3).move_to(DOWN*2.65)
    s.play(FadeIn(rows));s.wait(1)

    s.start('Across-chirp phase gives velocity — but it wraps',
        'Schematic candidate intersections, not measured ambiguity function. Uniform chirp starts separated by P give fslow modulo 1/P. Fast frequencies also wrap modulo fs. Approximate velocity alias interval c/(2fcP), about76m/s here. A fixed fast-frequency contour intersects each possible slow-frequency alias, leaving several range/velocity candidates. Bounds remove impossible candidates but do not guarantee uniqueness. A correction-bank full waveform comparison can reject approximate aliases when distinguishable residual phase exists, but cannot break exact sampled-data ambiguity. Additional PRI, slope diversity or external information is needed if multiple candidates remain equivalent.')
    axes=Axes(x_range=[.7,1.3,.2],y_range=[0,600,200],x_length=6.1,y_length=3.1,
        axis_config={'color':MUTED,'include_tip':False}).move_to([-2.7,.35,0])
    line=Line(axes.c2p(1.2,0),axes.c2p(.732,600),color=ORANGE,stroke_width=4)
    s.play(Create(axes),Create(line))
    s.play(FadeIn(s.prose('Range r0',22).next_to(axes,DOWN)),FadeIn(s.prose('Velocity v0',22).next_to(axes,UP)))
    levels=VGroup();dots=VGroup()
    for v in [40,116,192,268,344,420,496,572]:
        levels.add(DashedLine(axes.c2p(.7,v),axes.c2p(1.3,v),color=BLUE))
        dots.add(Dot(axes.c2p(1.2-.00078*v,v),color=PURPLE,radius=.065))
    s.play(LaggedStart(*[Create(a) for a in levels],lag_ratio=.15),run_time=2)
    s.play(FadeIn(dots))
    explanation=VGroup(s.prose('Blue: possible velocity aliases.',24,BLUE),
        s.prose('Purple: candidate solutions.',24,PURPLE),
        s.eq(r'f_{\rm slow}=\hat f_{\rm slow}+m/P',32),
        s.prose('m is an integer; P is chirp spacing.',22),
        s.eq(r'\Delta v\simeq\frac{c}{2f_cP}\approx76\ \mathrm{m/s}',30)).arrange(DOWN,buff=.32).move_to([3.6,.4,0])
    s.play(FadeIn(explanation))
    rows=VGroup(s.prose('Test every candidate inside the bounds against the full echo template.',25),
        s.prose('Bounds and waveform details can reject candidates; uniqueness is not guaranteed.',24),
        s.prose('Still ambiguous? Use different chirp spacings / slopes, or independent information.',24,BLUE)).arrange(DOWN,buff=.28).move_to(DOWN*2.65)
    s.play(FadeIn(rows));s.wait(1)

    s.start('Our choice: the sparse two-stage FFT',
        'Recommendation for this analyzer: keep the sparse two-stage implementation. Arithmetic model per RX and correction, not measured end-to-end runtime. K acquired chirps, Lf fast FFT length, Ls slow FFT length, M unique requested fast bins per correction group. Long comparator includes missing-time intervals and assumes an equivalent continuous-tone correction can be constructed. M329 is illustrative, not a measured group count. Eightfold padding gives K8,Lf2048,Ls64,Lgap32768: sparse306560 versus long491520 Llog2L units. Do not claim full-bank timing or extra SNR from factorization.')
    names=['Every chirp → FFT','Needed columns → FFT','Keep complex matched sums']
    boxes=VGroup(*[RoundedRectangle(width=3.9,height=1.1,stroke_color=BLUE).move_to([x,1.55,0]) for x in [-4.3,0,4.3]])
    s.play(FadeIn(boxes),*[FadeIn(s.prose(name,23,BLUE).move_to(box)) for name,box in zip(names,boxes)])
    s.play(*[GrowArrow(Arrow([x,1.55,0],[x+.35,1.55,0],buff=0,color=MUTED)) for x in [-2.3,2]])
    formula=s.eq(r'W=K L_f\log_2L_f+M L_s\log_2L_s',40).move_to(UP*.05)
    s.play(Write(formula))
    rows=VGroup(s.prose('K: chirps. Lf: fast FFT length. Ls: across-chirp FFT length.',24),
        s.prose('M: unique frequency columns requested by the physical grid.',24),
        s.prose('Example: 8 chirps, 329 columns → about 1.6 times less FFT work',26,BLUE),
        s.prose('than a gap-filled long FFT with comparable padding.',24,BLUE),
        s.prose('This is an arithmetic estimate; end-to-end runtime has not been compared.',23),
        s.prose('The saving is unused-column work and idle gaps, not a change in coherent SNR.',23)).arrange(DOWN,buff=.24).move_to(DOWN*1.95)
    s.play(FadeIn(rows));s.wait(1)
