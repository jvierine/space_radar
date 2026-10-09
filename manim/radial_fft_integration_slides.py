"""Full-vector versus separable FFT integration, by Juha Vierinen.
The cost comparison is analytical; no equivalent zero-filled FFT benchmark is claimed.
"""
from manim import *


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
    s.start('Two FFTs or one long FFT: which is faster?',
        'Costs per receiver per correction, omitting common preprocessing: fast FFTs K Lf log2 Lf, slow FFTs M Ls log2 Ls; M needed fast-frequency columns. One continuous-tone FFT with zeros at unobserved times costs Lgap log2 Lgap. These costs alone do not compare equivalent full range-velocity-acceleration banks: FMCW resets give different fast/slow frequencies, so a one-frequency approach needs additional range-dependent correction or another parameter loop. Long FFT uniform timing also needs special treatment since period*fs=317.125; ordinary zero insertion cannot encode fractional sample starts. No end-to-end comparative benchmark exists. Earlier timing compares one toy FFT with direct sums only. Both sums have same ideal SNR; temporal colored noise calls for C inverse q.')
    s.content([
        r'\text{two axes: }O(KL_f\log_2L_f+ML_s\log_2L_s)',
        r'\text{one gapped time axis: }O(L_{\rm gap}\log_2L_{\rm gap})',
    ],[
        'K: acquired chirps. Lf and Ls: padded FFT lengths. M: needed frequency columns.',
        'L(gap): padded length including the unobserved time intervals.',
        'Two axes avoid filling idle gaps and preserve the chirp-reset timing.',
        'One long FFT needs a continuous-tone model and a uniform time grid.',
        'Our chirp period is 317.125 sample intervals: simple zero insertion is insufficient.',
        'Which is faster for equivalent searches? Not benchmarked yet.',
        'The earlier toy timing compares FFT with direct sums, not these two methods.',
    ])
    s.start('An eight-chirp complexity estimate',
        'Illustrative arithmetic per RX, I/Q orientation and correction group. Actual code Nfast=nextpow2(8*225)=2048, Nslow=nextpow2(8*8)=64. Work units L log2 L: two-stage 180224+384M. At fs12.5MHz and P25.37us, gapped duration spans 7*317.125+225=2444.875 sample intervals. Eightfold padded long FFT nextpow2(8*2444.875)=32768, work491520. M329 is illustrative, not measured group size: cost306560, ratio1.603. All2048 columns cost966656, long cheaper1.967. Excludes scoring, corrections, fractional sample start handling and hardware overhead. Long FFT uses a different bank structure, so this is not an equivalent full-pipeline runtime estimate.')
    s.content([
        r'K=8,\quad N=225,\quad L_f=2048,\quad L_s=64',
        r'W_{\rm 2D}=180224+384M,\qquad W_{\rm long}=491520',
    ],[
        'Assume eightfold padding; count L log2(L) work, per RX and correction.',
        '329 needed columns: 306,560 units → about 1.6 times less FFT work.',
        'All 2,048 columns: 966,656 units → long FFT costs about half as much.',
        'Sparse second-stage columns are the useful computational saving.',
        'This estimates FFT arithmetic, not the full range / velocity / acceleration bank.',
        'Corrections, scoring, fractional timing and hardware overhead are excluded.',
    ])
