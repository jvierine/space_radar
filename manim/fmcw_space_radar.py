"""FYS-3000, 7 October 2026. Scientific products: space_radar_assets.py.
Source footer defaults on; SHOW_PROVENANCE=0 hides it.
"""
from pathlib import Path
import json
import h5py
import numpy as np
from PIL import Image
from manim import *
from planck_to_ktb import Text
import radar_equation_noise as base
from radar_equation_noise import RadarEquationNoise, BG,FG,MUTED,GOLD,SIGNAL,RECEIVE,POWER,RED
HERE=Path(__file__).resolve().parent
DATA=HERE/'assets/fmcw_lecture.h5'
CAPTURE=HERE.parent
GIF=HERE/'assets/core_dense_dx32_e_field_zoom.gif'
TI='https://www.ti.com/lit/ds/symlink/iwr1843.pdf'
BOARD='https://www.ti.com/lit/ug/spruim4b/spruim4b.pdf'
MIE='https://doc.comsol.com/6.3/doc/com.comsol.help.models.rf.rcs_sphere/rcs_sphere.html'
EXTRA={
'FMCW space radar':'FMCW transmits a frequency ramp and mixes the returned RF with the local transmit oscillator. Receiver output is downconverted complex I/Q. Frequency depends on both delay and Doppler; a rapid trajectory needs the complete phase. Source: '+TI,
'Fraunhofer test':'Schematic, not surveyed geometry. Radar1 is at the left end, looking right toward the impact target, following the user location correction. Recorded shot 6606: 5 mm, 6700 m/s, 78 GHz, nominally toward impact target; radar1 offset 11–12 cm from projectile path and 9 cm before back wall. Source angle 82.7 degrees has unspecified reference, omitted. Shot 6607: 6 mm, 6285 m/s, 79.5 GHz, perpendicular setup, separate day; transverse distance not surveyed. Chamber 100 mbar. Source: capture NetCDF comments and IDs in HDF5. One metre drawn as assumed observed trajectory window, not vessel length.',
'Space radar parameters':'128 complex samples at 12.5 Msps = 10.24 us live ADC time. Chirp period 13.71 us inferred from reader timing, not measured timestamp spacing. Frame gaps omitted. Maximum 12 dBm per TX confirmed by user. 6/10 dBi gains are sensitivity assumptions. Board peak >10.5 dBi does not determine gain toward a particular target. Sources: '+TI+' ; '+BOARD,
'Inverse square spreading':'Spherical surface area 4*pi*R^2. Fixed directional gain is assumed. Far-field propagation. Source: '+base.RADAR_SOURCE,
'Two propagation legs':'Outgoing flux proportional R^-2; equivalent backscatter flux R^-2; receiving aperture fixed gives R^-4. Same target, wavelength, gains and losses. Source: '+base.RADAR_SOURCE,
'Why proximity matters':'Ratio (1000 km/1 km)^4=10^12=120 dB. Far field and unchanged target/aspect/gains assumed.',
'Antenna beam':'Gain=efficiency*directivity. Dmax=4*pi/Omega_A where Omega_A integrates normalized antenna power pattern over whole sphere. Uniform circular cone shown only as an idealized lossless pattern. Actual xWR1843BOOST fan-shaped pattern about +/-60 deg azimuth and +/-15 deg elevation. FoV definition differs from directivity solid angle. Sources: '+base.NOISE_SOURCE+' ; '+BOARD,
'Noise figure to temperature':'F=10^(NF/10), Te=(F-1)*290 K, Tsys=Tant+Te. At Tant=290 K and NF15 dB, Te8880.6 K, Tsys9170.6 K. We round to10000 K. NF depends on RXgain. Noise reference plane is antenna input. Source: '+TI,
'Rayleigh and Mie scattering':'Perfectly conducting sphere. Diameter/wavelength, not radius/wavelength. Rayleigh sigma=9*pi*a^2*(ka)^4. Exact Mie valid for all sizes. Rayleigh breakdown near d/lambda~0.1 is approximate. Source: '+MIE+' ; capture_analysis/analyze_captures.py sphere().',
'Millimetre wavelengths':'3 mm PEC sphere RCS from exact Mie solver. Comparing RCS alone; a complete link also contains lambda squared and antenna gains. 78 GHz sigma1.30024e-5 m^2 vs3GHz5.02513e-9, factor2587.5.',
'FDTD scattering demonstration':'Existing sabmod animation: meteor plasma at53.5 MHz, dense-core visualization, not simulation of the metal pellet or the IWR1843 band. Actual total/scattered Ez; coarse dx0.03125 m; unresolved brightest core not a calibrated RCS result. Scripts animate_e_field.py, maarsy_rcs.py, fdtd/src/main.rs; README_rcs.md smaller-box denser-plasma example. Source GIF '+str(GIF),
'Three simulated chirps':'Exact memo Eq25 phase model, three consecutive first chirps of illustrative parallel trajectory. Signal uses retarded range, quadratic delay term and per-chirp fast-time reset, HPF attenuation and assumed IF sideband. Input-equivalent sqrt(W), not calibrated ADC volts. Geometry not surveyed. Source: '+str(CAPTURE/'capture_analysis/analyze_captures.py'),
'Analogue chirp downconversion':'The receiver multiplies the echo by the conjugate local transmit chirp and low-pass filters the RF sum-frequency products. The residual complex beat signal contains range-delay and Doppler contributions; it is not the 78 GHz carrier. I and Q are its real and imaginary components. Ideal mixer phase is Eq25. The familiar fIF=-2SR/c+fD and fD=-2vr/lambda are approximations; the plotted fast moving-target waveform uses full retarded phase and the assumed receiver response. ADC samples are separated by idle gaps. Source: '+TI+' ; '+str(CAPTURE/'capture_analysis/analyze_captures.py'),
'Trying matched filters':'Candidates change speed while preserving initial geometry. Exact complex IQ templates include amplitude and phase. Correct filter product s*q_conj is real nonnegative; wrong product oscillates and cancels. Integral uses actual samples at fs12.5MHz. Normalized squared coherence shown, not measured detection SNR.',
'Energy and bandwidth':'Optimal matched-filter white-noise SNR Es/(kB*Tsys). T_int denotes the summed duration of samples used coherently, excluding ADC gaps. Equivalent constant-amplitude uniform estimator has two-sided noise-equivalent bandwidth B=1/T_int. The full elapsed coherent window includes gaps and differs from T_int; phase must remain matched across that window. For variable amplitude use matched weights. Gaps do not collect echo energy.',
'Interactive simulator':'Open https://juha.no/fmcw/sim/ using the clickable simulator link in the browser deck. Explore target trajectory, 1–10 chirps including ADC padding, complex beat voltage, instantaneous frequency and spectrogram. The two measured template-search diagnostic slides were removed at the user request; underlying analysis products are retained.',
'One metre field of view':'Assumed useful trajectory length1m, not measured beam footprint. Tvis=1m/v; Nperiods=Tvis/13.71us. Shot6606:149.25us,10.8865periods(~11chirps);6607:159.11us,11.6053periods(~12chirps). Entry/exit chirps can be partial; counts are approximate period-equivalent dwell, not guaranteed full-chirp counts. Ten complete recorded chirps span133.63us and collect102.4usADC, so a10chirp window fits within the assumed dwell at BOTH recorded speeds. At7km/s:142.86us,10.42periods; at14km/s:71.43us,5.21periods. Actual visible path depends on beam, range, pointing, trajectory and sensitivity. Timing inferred from reader.',
'Three millimetre SNR':'Ideal thermal-limited, PEC3mm sphere,78GHz,12dBm,oneTXoneRX,gains6dBi each,Tsys10000K,loss0dB,constant range and correct phase. N chirps collect N*10.24us. Numbers are not calibrated measured SNR. Gains10dBi each add8dB.',
'Fraunhofer SNR estimates':'Conditional thermal budget at fixed1m, using actual projectile diameters and recorded frequencies. Maximum12dBm user-confirmed; gains6dBi each andTsys10000K assumed. Source HDF fmcw_lecture.h5. Actual geometry, pointing, varying range, HPFs, sideband selection, losses and clutter can change result.',
'Measured undecoded voltage':'Shot6607, RX0, frames34–37. Real part of stored complex IQ in ADC counts. No decoding or voltage calibration. Stable early fast-time structure interpreted as stationary vessel walls/objects (and possible fixed leakage). After disturbance rapidly changing strong background. Oscillating chamber is the physical interpretation proposed by user; data alone do not prove origin. Time axis uses assumed13.71us, frame gaps omitted.',
'Quiet background subtraction':'Complex mean per fast-time bin over chirps0–1623 only, then subtract the same mean from every chirp. No matched filter/FFT/decoding in this real-component display. Product: output/background_subtracted_6607_frames_34_35_36_37.h5; source capture_analysis/plot_background_subtracted_views.py. Changing background remains after fixed mean removal.',
'Disturbance and self noise':'Time-varying chamber echoes are structured clutter, called self noise here. They can dominate a weak pellet echo; thermal SNR budget alone does not describe practical detectability. Real component amplitude image is not itself a calibrated echo-power measurement. Use the pre-disturbance interval and phase-aware templates; do not claim pellet detection from this display.'
}
base.NOTES.update(EXTRA)
base.NOTES.update({
    'Trajectory coordinates':'Radar at(x0,y0), projectile at(v0*t,0); time origin is first chirp of EACH fitted window. x0 along track, y0 perpendicular; rangeR=sqrt((v0*t-x0)^2+y0^2). Signed v0 is along track, not instantaneous radial speed. Reversing(x0,v0) leaves range history unchanged in this isotropic monostatic model.'
    ,'Measured frequency interpretation':'Broad signed-velocity bank(-7000..7000m/s,100m/s grid and local refinement), all strict quiet windows. Carrier Doppler=-f0*2vr/(c+vr); range beat=-S*tau; total IF includes chirp-motion term too. Strict frame34 maximum at1.098km/s is a background template fit near longitudinal boundary, not an identified object. Physical sign and stored IQ conjugation are separate. Sources search_signed_velocities.py, show_signed_velocity.py.'
})
base.NOTES.update({
    'Aperture numbers':'Ae=G lambda squared/(4pi), at78GHz, gain6dBi gives4.68mm2, gain10dBi11.76mm2. These are effective collecting areas, not board physical areas. Efficiency contained in gain. Source: MIT6.013 antenna reciprocity.',
    'Accumulated filter output':'Plot absolute cumulative complex integral signal*conjugate template, normalized by full signal energy. Actual Eq25 samples. Correct template yields nonnegative products and full energy at end; wrong templates largely cancel. Trial6200/7200/6700m/s. Normalization does not mean the accumulated quantity is instantaneous power.'
    ,'Pre-disturbance measured search':'Actual shot6607 complex-mean-subtracted RX0. Signed velocities-7..7km/s in100m/ssteps, localrefinement; x0-1.2..1.2m,y0.05..1m. All140643overlapping8chirpwindows in frames0–33and34before1624; bothIQsigns; no whitening or frame stitching. Negative bank exact mirror(x0,v0)->(-x0,-v0). Strict change audit first sustained power increase median+6MADscale for3chirps. Same plotted mean retained. Strictframe34maximum chirp289, score71.79, x0=1.2m(longitudinal boundary),v0=1.0976km/s,y0=.6468m, explains5.13%; stationary templates dominate earlycontrols. Not an isolated fast echo. Score scalarresidualvariance-normalized, notthermalSNR/falsealarmprobability. Source search_signed_velocities.py,show_signed_velocity.py; CPU/MPS verified, broadinjectioncoherence.9947.'
})

class FMCWSpaceRadar(RadarEquationNoise):
    def setup(self):
        super().setup()
        self.provenance.become(Text('Sources: fmcw_space_radar.py; space_radar_assets.py',font_size=13,color=MUTED).to_corner(DR,buff=.12))
        self.h=h5py.File(DATA,'r')

    def construct(self):
        self.introduction();self.test_setup();self.parameters();self.coordinates();self.spreading();self.two_legs();self.proximity()
        self.radar_equation();self.transmit();self.cross_section();self.beam();self.aperture();self.aperture_numbers()
        self.scattering();self.wavelength();self.fdtd();self.noise();self.temperature()
        self.downconversion();self.three_chirps();self.filters();self.accumulated();self.energy_bandwidth();self.filter_response();self.bandwidth()
        self.integration_time();self.visibility();self.snr_table();self.fraunhofer_snr()
        self.measured();self.subtracted()
        for shot,frames in [('6607',[34,35,36,37]),('6606',[43,44,45,46])]:
            self.chirp_spectra(shot,frames,False);self.chirp_spectra(shot,frames,True)
        self.clutter();self.matched_noise_definition();self.four_chirp_search();self.matched_noise_controls();self.simulator()

    def text_slide(self,name,title,lines,formula=None):
        self.start_slide(name);self.play(FadeIn(self.title(title)))
        y=1.7
        if formula:
            self.play(Write(self.equation(formula,size=48,color=GOLD).move_to(UP*1.8)));y=.5
        for line in lines:
            self.play(FadeIn(self.prose(line,29).move_to(UP*y)),run_time=.5);y-=.85
        self.finish()

    def introduction(self):
        self.start_slide('FMCW space radar');self.play(FadeIn(self.title('FMCW space radar: transmit, scatter, mix')))
        geometry=self.geometry().shift(UP*.8);self.play(FadeIn(geometry))
        pulse=Dot(geometry[2].get_start(),color=POWER,radius=.13)
        self.add(pulse);self.play(MoveAlongPath(pulse,geometry[2]),run_time=1.1,rate_func=linear)
        self.play(pulse.animate.set_color(RECEIVE).move_to(geometry[3].get_start()),run_time=.15)
        self.play(MoveAlongPath(pulse,geometry[3]),run_time=1.1,rate_func=linear);self.remove(pulse)
        ax=Axes(x_range=[0,3,1],y_range=[0,3,1],x_length=5,y_length=1.7,tips=False,axis_config={'color':MUTED}).move_to([-3,-1.2,0])
        tx=VGroup(*[ax.plot_line_graph([i,i+.9],[.2,2.7],add_vertex_dots=False,line_color=POWER) for i in range(3)])
        rx=VGroup(*[ax.plot_line_graph([i+.15,i+1.05],[.2,2.7],add_vertex_dots=False,line_color=RECEIVE) for i in range(3)])
        self.play(Create(ax),Create(tx));self.play(Create(rx))
        self.play(FadeIn(self.prose('Frequency ramps',25).next_to(ax,DOWN)),Write(self.equation(r'z(t)=\mathrm{LPF}\{s_{\rm echo}(t)s_{\rm TX}^{*}(t)\}',size=31).move_to([3.3,-.9,0])))
        self.play(FadeIn(self.prose('Delay and motion become phase and frequency in complex I/Q',29).move_to(DOWN*2.85)))
        self.finish()

    def test_setup(self):
        self.start_slide('Fraunhofer test');self.play(FadeIn(self.title('Fraunhofer hypervelocity chamber tests')))
        wall=RoundedRectangle(width=10.7,height=2.5,corner_radius=.18,color=MUTED).move_to(UP*.7)
        path=Arrow([-5,.7,0],[4.7,.7,0],color=GOLD,buff=0)
        impact=Line([4.8,-.4,0],[4.8,1.8,0],stroke_width=9,color=RED)
        r1=Triangle(color=SIGNAL,fill_opacity=.5).scale(.18).rotate(-PI/2).move_to([-4.9,-.15,0])
        r2=Triangle(color=RECEIVE,fill_opacity=.5).scale(.18).move_to([.3,-.85,0])
        self.play(Create(wall),Create(path),Create(impact),FadeIn(r1),FadeIn(r2))
        pellet=Dot([-4.8,.7,0],radius=.085,color=GOLD);self.add(pellet)
        self.play(pellet.animate.move_to([4.7,.7,0]),run_time=1.3,rate_func=linear)
        labs=VGroup(self.prose('Projectile',23,GOLD).move_to([-3,1.4,0]),self.prose('Impact target',22,RED).move_to([4.8,2.3,0]),self.prose('Radar 1: toward target',22,SIGNAL).move_to([-3.7,-.65,0]),self.prose('Radar 2: perpendicular',22,RECEIVE).move_to([-.8,-1.5,0]))
        self.play(FadeIn(labs))
        for i,line in enumerate(['6606: 5 mm pellet, 6.70 km/s; 6607: 6 mm pellet, 6.285 km/s',
                                 'Radar 1 at the left end, looking along the projectile path',
                                 '100 mbar; separate shots on separate days; schematic, not to scale']):
            self.play(FadeIn(self.prose(line,25,MUTED if i==2 else FG).move_to(DOWN*(2.05+.52*i))))
        self.finish()

    def parameters(self):
        self.text_slide('Space radar parameters','What was recorded',[
            'IWR1843: 78 / 79.5 GHz; wavelength about 3.8 mm',
            'Chirp slope: 100 / 60 MHz per microsecond',
            '128 complex ADC samples at 12.5 MS/s: 10.24 microseconds',
            'Chirp period about 13.71 microseconds; 4096 chirps per frame',
            'Max transmit power: 12 dBm = 15.85 mW per TX',
            'Budget: one TX, one RX; 6 and 10 dBi directional gain assumptions'])

    def coordinates(self):
        self.start_slide('Trajectory coordinates');self.play(FadeIn(self.title('Position along the path is different from radar range')))
        path=Arrow([-5,-.8,0],[5.5,-.8,0],buff=0,color=MUTED)
        pellet=Dot([-3,-.8,0],color=GOLD,radius=.12);radar=Dot([1,1.5,0],color=SIGNAL,radius=.12)
        longitudinal=DoubleArrow([-3,-1.4,0],[1,-1.4,0],buff=0,color=GOLD)
        transverse=DoubleArrow([1,-.8,0],[1,1.5,0],buff=.1,color=SIGNAL)
        slant=Line([-3,-.8,0],[1,1.5,0],color=RECEIVE)
        velocity=Arrow([-3,-.5,0],[-1.1,-.5,0],buff=0,color=POWER)
        self.play(Create(path),FadeIn(pellet),FadeIn(radar),Create(slant),Create(longitudinal),Create(transverse),Create(velocity))
        labs=VGroup(self.prose('Projectile at window start',24,GOLD).move_to([-3.5,-2,0]),self.prose('Radar',25,SIGNAL).next_to(radar,UP),
                    self.equation(r'x_0',size=36,color=GOLD).move_to([-1,-1.1,0]),self.equation(r'y_0',size=36,color=SIGNAL).move_to([1.5,.4,0]),
                    self.equation(r'R(0)',size=36,color=RECEIVE).move_to([-1.5,.65,0]),self.equation(r'v_0',size=34,color=POWER).move_to([-2,-.05,0]))
        self.play(FadeIn(labs))
        self.play(Write(self.equation(r'R(t)=\sqrt{(v_0t-x_0)^2+y_0^2},\qquad R(0)=\sqrt{x_0^2+y_0^2}',size=40).move_to(DOWN*2.75)))
        self.play(FadeIn(self.prose('x0 along the path; y0 perpendicular; v0 along the path, not radial',23,MUTED).move_to(DOWN*3.4)))
        self.finish()

    def spreading(self):
        self.start_slide('Inverse square spreading');self.play(FadeIn(self.title('The same power spreads over a larger surface')))
        centre=np.array([-3,0,0]);origin=Dot(centre,color=POWER)
        rings=VGroup(*[Circle(radius=r,color=POWER,stroke_opacity=.8).move_to(centre) for r in [.65,1.3,2.6]])
        self.add(origin);self.play(LaggedStart(*[Create(r) for r in rings],lag_ratio=.3))
        self.play(Write(self.equation(r'A_{\rm sphere}=4\pi R^2',size=45).move_to([3,1.1,0])))
        self.play(Write(self.equation(r'S_{\rm inc}=\frac{P_tG_t}{4\pi R^2}',size=49,color=POWER).move_to([3,-.3,0])))
        self.play(FadeIn(self.prose('Double the range: four times the area, one quarter the flux',29).move_to(DOWN*3)))
        self.finish()

    def two_legs(self):
        self.start_slide('Two propagation legs');self.play(FadeIn(self.title('The echo spreads again on its way back')))
        self.play(FadeIn(self.geometry().shift(UP*1.2)))
        self.play(Write(self.equation(r'S_{\rm inc}=\frac{P_tG_t}{4\pi R^2}',size=44,color=POWER).move_to([-3,-.55,0])))
        self.play(Write(self.equation(r'S_{\rm back}=S_{\rm inc}\frac{\sigma}{4\pi R^2}',size=44,color=RECEIVE).move_to([3,-.55,0])))
        self.play(Write(self.equation(r'P_r=S_{\rm back}A_e/L\quad\propto\quad\frac{1}{R^2}\frac{1}{R^2}=\frac{1}{R^4}',size=48,color=GOLD).move_to(DOWN*1.65)))
        self.play(FadeIn(self.prose('Illumination × scattering × receiving aperture',29).move_to(DOWN*2.9)))
        self.finish()

    def proximity(self):
        self.text_slide('Why proximity matters','Get close to a small target',[
            'The same target at 1 km gives a trillion times more echo power',
            'than at 1000 km, with the same antennas and transmitted power.',
            'This is a 120 dB advantage in received power.'],r'\frac{P_r(1\,\mathrm{km})}{P_r(1000\,\mathrm{km})}=1000^4=10^{12}')

    def beam(self):
        self.start_slide('Antenna beam');self.play(FadeIn(self.title('Antenna gain concentrates power into a beam')))
        for x,gain,angle,col in [(-3,6,120,SIGNAL),(3,10,73.74,RECEIVE)]:
            half=angle/2*DEGREES;apex=np.array([x,-.5,0])
            points=[apex,apex+2.7*np.array([np.sin(half),np.cos(half),0]),apex+2.7*np.array([-np.sin(half),np.cos(half),0])]
            cone=Polygon(*points,color=col,fill_opacity=.18)
            self.play(Create(cone),FadeIn(self.prose(f'{gain} dBi: ideal cone {angle:.0f} degrees',25,col).move_to([x,2.2,0])))
        self.play(Write(self.equation(r'G=\eta D,\qquad D_{\max}=\frac{4\pi}{\Omega_A},\qquad \Omega_A=\int p_n(\theta,\varphi)\,d\Omega',size=37).move_to(DOWN*1.55)))
        self.play(FadeIn(self.prose('Ideal lossless cones; the board has a fan-shaped beam',25,MUTED).move_to(DOWN*2.45)))
        self.play(FadeIn(self.prose('Board: about ±60 degrees azimuth, ±15 degrees elevation; peak >10.5 dBi',23,MUTED).move_to(DOWN*3.05)))
        self.finish()

    def scattering(self):
        self.start_slide('Rayleigh and Mie scattering');self.play(FadeIn(self.title('Small spheres: Rayleigh to Mie scattering')))
        g=self.h['scattering'];x=np.log10(g['diameter_over_wavelength'][:]);mie=np.log10(g['mie_normalized_rcs'][:]);ray=np.log10(g['rayleigh_normalized_rcs'][:])
        ax=Axes(x_range=[-3,1,1],y_range=[-9,2,2],x_length=10.0,y_length=4,tips=False,axis_config={'color':MUTED}).move_to(RIGHT*.4+UP*.05)
        self.play(Create(ax))
        for y,col in [(ray,POWER),(mie,SIGNAL)]:
            mask=y<=2;curve=VMobject(color=col).set_points_as_corners([ax.c2p(a,b) for a,b in zip(x[mask],y[mask])]);self.play(Create(curve),run_time=1)
        for xx,lab in [(-3,'0.001'),(-2,'0.01'),(-1,'0.1'),(0,'1'),(1,'10')]:self.add(self.prose(lab,19,MUTED).move_to(ax.c2p(xx,-9)+DOWN*.22))
        y_ticks=VGroup(*[self.equation(r'10^{'+str(yy)+'}',size=20,color=MUTED).next_to(ax.c2p(-3,yy),LEFT,buff=.15) for yy in [-8,-6,-4,-2,0,2]])
        y_label=self.prose('RCS / geometric area (log scale)',23,MUTED).rotate(PI/2).next_to(y_ticks,LEFT,buff=.30)
        self.add(y_ticks,y_label,self.prose('Diameter / wavelength (log scale)',24).next_to(ax,DOWN,buff=.5))
        self.play(FadeIn(self.prose('Rayleigh approximation',23,POWER).move_to([-3,2.6,0])),FadeIn(self.prose('Exact conducting-sphere Mie',23,SIGNAL).move_to([3,2.6,0])))
        self.play(FadeIn(self.prose('Around diameter = 0.1 wavelength, use the full scattering solution',25).move_to(DOWN*3.1)))
        self.finish()

    def aperture_numbers(self):
        ae=self.h['antenna/effective_aperture_mm2'][:]
        self.text_slide('Aperture numbers','From gain to collecting area',[
            'At 78 GHz the wavelength is 3.84 mm.',
            f'6 dBi gain: G = 3.98, effective aperture = {ae[0]:.2f} square millimetres',
            f'10 dBi gain: G = 10, effective aperture = {ae[1]:.2f} square millimetres',
            'Effective aperture describes collected power; it is not the board area.',
            'Efficiency is already in gain; additional loss is represented by L.'],r'A_e=\frac{G_r\lambda^2}{4\pi},\qquad P_r=S_{\rm back}A_e/L')

    def wavelength(self):
        sigma=self.h['scattering/three_mm_rcs_m2'][:]
        self.text_slide('Millimetre wavelengths','Why millimetre wavelengths for millimetre targets?',[
            f'3 mm conducting sphere at 3 GHz: RCS = {sigma[0]:.2g} square metres',
            f'At 78 GHz: RCS = {sigma[2]:.2g} square metres',
            f'The backscatter cross section is about {sigma[2]/sigma[0]:.0f} times larger.',
            'The full radar budget also includes wavelength squared and antenna gain.'],r'\sigma_{\rm Rayleigh}=9\pi a^2(ka)^4,\qquad k=2\pi/\lambda')

    def fdtd(self):
        self.start_slide('FDTD scattering demonstration');self.play(FadeIn(self.title('Watch an incident wave scatter')))
        im=Image.open(GIF);im.seek(0)
        mob=ImageMobject(im.convert('RGB')).scale_to_fit_height(5).move_to(UP*.1);self.add(mob)
        for j in range(0,im.n_frames,5):
            im.seek(j);mob.become(ImageMobject(im.convert('RGB')).scale_to_fit_height(5).move_to(UP*.1));self.wait(.2)
        self.play(FadeIn(self.prose('Existing FDTD meteor-plasma example at 53.5 MHz; illustrates scattering',25).move_to(DOWN*2.8)))
        self.play(FadeIn(self.prose('Different target and wavelength from the chamber pellet',23,MUTED).move_to(DOWN*3.3)))
        self.finish()

    def temperature(self):
        self.text_slide('Noise figure to temperature','From receiver noise figure to system temperature',[
            'NF = 15 dB gives F = 31.62, using reference temperature 290 K.',
            'Receiver equivalent noise temperature: 8881 K',
            'With a 290 K antenna: system temperature = 9171 K',
            'We use 10 000 K as a rounded budget assumption.'],r'F=10^{\mathrm{NF}/10},\quad T_e=(F-1)T_0,\quad T_{\rm sys}=T_{\rm ant}+T_e')

    def sampled_curve(self,ax,t,values,col):
        return VGroup(*[VMobject(color=col,stroke_width=2).set_points_as_corners([ax.c2p(a,b) for a,b in zip(t[k],values[k])]) for k in range(len(t))])

    def downconversion(self):
        self.start_slide('Analogue chirp downconversion');self.play(FadeIn(self.title('The receiver measures a downconverted beat signal')))
        self.play(FadeIn(self.prose('Mix the delayed echo with the conjugate of the local transmit chirp.',27).move_to(UP*2.25)))
        self.play(Write(self.equation(r'z_{\rm IF}(t)=\mathrm{LPF}\{s_{\rm echo}(t)s_{\rm TX}^{*}(t)\}=I(t)+iQ(t)',size=37).move_to(UP*1.35)))
        self.play(FadeIn(self.prose('The RF carrier is removed; range delay and motion leave oscillations.',27).move_to(UP*.45)))
        self.play(Write(self.equation(r'z_{\rm IF}=A(t)e^{i\phi(t)},\qquad f_{\rm IF}=\frac{1}{2\pi}\frac{d\phi}{dt}',size=38).move_to(DOWN*.5)))
        self.play(Write(self.equation(r'f_{\rm IF}\simeq-\frac{2SR}{c}+f_D,\qquad f_D\simeq-\frac{2v_r}{\lambda}',size=37).move_to(DOWN*1.55)))
        self.play(FadeIn(self.prose('Range beat and carrier Doppler both contribute, in this I/Q convention.',25,MUTED).move_to(DOWN*2.45)))
        self.play(FadeIn(self.prose('For a rapidly moving target, use the full Eq. 25 phase on the next slide.',25,GOLD).move_to(DOWN*3.1)))
        self.finish()

    def three_chirps(self):
        self.start_slide('Three simulated chirps');self.play(FadeIn(self.title('Three chirps of the downconverted beat signal')))
        g=self.h['matched_filter_demo'];t=g['time_s'][:]*1e6;s=g['signal_iq'][:];s=s/abs(s).max()
        self.play(Write(self.equation(r'\phi=2\pi[-f_0\tau-Su\tau+\tfrac12S\tau^2],\quad\tau=2R_{\rm ret}/c',size=33).move_to(UP*2.35)))
        self.play(FadeIn(self.prose('Residual echo phase; u resets each chirp; tau is the round-trip delay.',23,MUTED).move_to(UP*1.7)))
        ax=Axes(x_range=[0,38,10],y_range=[-1,1,1],x_length=10.5,y_length=2.3,tips=False,axis_config={'color':MUTED}).move_to(DOWN*.45)
        self.play(Create(ax),Create(self.sampled_curve(ax,t,s.real,SIGNAL)),Create(self.sampled_curve(ax,t,s.imag,RECEIVE)),run_time=1.5)
        self.add(self.prose('I: real component',23,SIGNAL).move_to([-3,1,0]),self.prose('Q: imaginary component',23,RECEIVE).move_to([3,1,0]))
        for tick in [0,10,20,30]:self.add(self.prose(str(tick),20,MUTED).move_to(ax.c2p(tick,-1)+DOWN*.22))
        y_ticks=VGroup(*[self.prose(str(tick),19,MUTED).next_to(ax.c2p(0,tick),LEFT,buff=.12) for tick in [-1,0,1]])
        self.add(y_ticks,self.prose('Normalized I, Q',21,MUTED).rotate(PI/2).next_to(y_ticks,LEFT,buff=.2))
        self.play(FadeIn(self.prose('Elapsed time (microseconds)',23).next_to(ax,DOWN,buff=.5)))
        self.play(FadeIn(self.prose('128 complex samples per chirp; 10.24 microseconds ADC; 13.71 microseconds period',23,MUTED).move_to(DOWN*2.65)))
        self.play(FadeIn(self.prose('30.72 microseconds of data across 37.66 microseconds; no samples in the gaps',23,MUTED).move_to(DOWN*3.2)))
        self.finish()

    def filters(self):
        self.start_slide('Trying matched filters');self.play(FadeIn(self.title('Try a phase history: do the oscillations cancel?')))
        g=self.h['matched_filter_demo'];t=g['time_s'][:]*1e6;mixed=g['mixed_iq'][:];ints=g['cumulative_integral'][:];v=g['candidate_velocity_m_s'][:];coh=g['power_coherence'][:]
        ax=Axes(x_range=[0,38,10],y_range=[-1,1,1],x_length=6.8,y_length=3,tips=False,axis_config={'color':MUTED}).move_to([-2.4,0,0])
        self.add(ax,self.prose('Real part of signal × conjugate template',22).next_to(ax,DOWN))
        current=None;label=None;result=None
        for j,col in enumerate([RED,POWER,RECEIVE]):
            curve=self.sampled_curve(ax,t,mixed[j].real/abs(mixed).max(),col)
            lab=self.prose(f'Trial speed: {v[j]/1000:.2f} km/s',29,col).move_to(UP*2.3)
            res=VGroup(self.prose('Normalized match power',24),self.equation(f'{coh[j]:.4f}',size=52,color=col),self.prose('Correct phase adds coherently' if j==2 else 'Oscillations cancel in the integral',24,col)).arrange(DOWN,buff=.5).move_to([4,0,0])
            if current:self.play(FadeOut(current),FadeOut(label),FadeOut(result))
            self.play(Create(curve),FadeIn(lab),FadeIn(res),run_time=1);self.wait(1)
            current,label,result=curve,lab,res
        self.play(Write(self.equation(r'M(\theta)=\int z(t)q_\theta^*(t)\,dt',size=38).move_to(DOWN*2.8)))
        self.finish()

    def energy_bandwidth(self):
        self.text_slide('Energy and bandwidth','Two equivalent views of coherent sensitivity',[
            'Energy view: match the waveform and collect the echo energy.',
            'Filter view: preserve the aligned signal and narrow the noise bandwidth.',
            'For constant amplitude, echo energy = received power × integration time.',
            'Motion and changing amplitude require the corresponding matched weights.'],r'\rho=\frac{E_s}{k_BT_{\rm sys}}=\frac{P_rT_{\rm int}}{k_BT_{\rm sys}}=\frac{P_r}{k_BT_{\rm sys}B_{\rm eff}}')

    def accumulated(self):
        self.start_slide('Accumulated filter output');self.play(FadeIn(self.title('The correct filter accumulates echo energy')))
        g=self.h['matched_filter_demo'];t=g['time_s'][:].ravel()*1e6
        integral=g['cumulative_integral'][:];energy=g.attrs['signal_energy_J']
        ax=Axes(x_range=[0,38,10],y_range=[0,1.05,.5],x_length=10.7,y_length=3.5,tips=False,
                axis_config={'color':MUTED},x_axis_config={'include_numbers':True},y_axis_config={'include_numbers':True}).move_to(DOWN*.1)
        self.play(Create(ax))
        for j,col in enumerate([RED,POWER,RECEIVE]):
            curve=VMobject(color=col,stroke_width=3).set_points_as_corners([ax.c2p(x,y) for x,y in zip(t,abs(integral[j])/energy)])
            self.play(Create(curve),FadeIn(self.prose(f'{g["candidate_velocity_m_s"][j]/1000:.2f} km/s',25,col).move_to([-3.5+j*3.5,2.25,0])),run_time=1)
        self.play(FadeIn(self.prose('Time (microseconds)',25).next_to(ax,DOWN,buff=.16)),Write(self.equation(r'\left|\int_0^t z(t\prime)q_\theta^*(t\prime)\,dt\prime\right|/E_s',size=31).move_to(DOWN*3.1)))
        self.finish()

    def integration_time(self):
        self.text_slide('Ten seconds of coherent measurements','The bandwidth of your coherent filter',[
            'Ten seconds of continuous, phase-matched data: bandwidth 0.1 Hz.',
            'Integration time counts the samples used coherently, excluding gaps.',
            'Three chirps: integration time 30.72 microseconds; span 37.66 microseconds.',
            'Their constant-amplitude averaging bandwidth is 32.55 kHz.',
            'Maintain the matched phase across the entire window, including gaps.'],r'B_{\rm eff}=1/T_{\rm int},\qquad T_{\rm int}=N_{\rm samples}/f_s')

    def visibility(self):
        self.start_slide('One metre field of view');self.play(FadeIn(self.title('Expect about 11--12 chirps from the pellet')))
        window=Rectangle(width=9,height=.65,color=RECEIVE,fill_opacity=.1).move_to(UP*1.95)
        self.play(Create(window),FadeIn(self.prose('Assumed 1 metre of visible trajectory',26).move_to(UP*2.65)))
        dot=Dot([-4.5,1.95,0],color=GOLD);self.add(dot)
        self.play(dot.animate.move_to([4.5,1.95,0]),run_time=1.5,rate_func=linear)
        self.play(Write(self.equation(r'T_{\rm vis}=\ell/v,\qquad N_{\rm chirps}\simeq\frac{T_{\rm vis}}{T_{\rm rep}},\qquad T_{\rm rep}=13.71\ \mu\mathrm s',size=35).move_to(UP*1.05)))
        cells=[[self.prose(x,26,GOLD if i==0 else FG) for x in row] for i,row in enumerate([
            ['Shot / speed','1 m transit','Chirp periods'],
            ['6606 / 6.70 km/s','149.3 microseconds','10.9 (about 11)'],
            ['6607 / 6.285 km/s','159.1 microseconds','11.6 (about 12)']])]
        table=MobjectTable(cells,include_outer_lines=True,h_buff=.55,v_buff=.18,line_config={'color':MUTED,'stroke_width':1})
        table.scale_to_fit_width(10.8);table.scale_to_fit_height(1.75);table.move_to(DOWN*.35)
        self.play(Create(table),run_time=1.1)
        self.play(FadeIn(self.prose('A 10-chirp window fits: 133.63 microseconds span, 102.4 microseconds live ADC',24,RECEIVE).move_to(DOWN*1.7)))
        self.play(FadeIn(self.prose('At 7 km/s: about 10 chirps; at 14 km/s: about 5 chirps over 1 metre',24).move_to(DOWN*2.35)))
        self.play(FadeIn(self.prose('Entry/exit chirps can be partial; actual visibility depends on the beam and range',22,MUTED).move_to(DOWN*3)))
        self.finish()

    def number_table(self,headers,rows,y=.1):
        cells=[[self.prose(str(x),27,GOLD if i==0 else FG) for x in row] for i,row in enumerate([headers]+rows)]
        table=MobjectTable(cells,include_outer_lines=True,h_buff=.65,v_buff=.3,line_config={'color':MUTED,'stroke_width':1})
        table.scale_to_fit_width(10.8)
        if table.height>3.6:table.scale_to_fit_height(3.6)
        table.move_to(UP*y);self.play(Create(table),run_time=1.3)

    def snr_table(self):
        self.start_slide('Three millimetre SNR');self.play(FadeIn(self.title('3 mm sphere: ideal thermal SNR (dB)')))
        g=self.h['three_mm'];a=g['snr_gain_6_db'][:];n=g['chirps'][:]
        self.number_table(['Chirps','1 m','10 m','100 m'],[[str(k)]+[f'{v:.1f}' for v in row] for k,row in zip(n,a)],.25)
        self.play(FadeIn(self.prose('78 GHz; 12 dBm TX; 6 dBi TX and RX; 10 000 K; one RX',25).move_to(DOWN*2.1)))
        self.play(FadeIn(self.prose('10 dBi on both antennas adds 8 dB to every entry',27,GOLD).move_to(DOWN*2.75)))
        self.play(FadeIn(self.prose('Constant range, correct phase, conducting sphere; no added losses or clutter',22,MUTED).move_to(DOWN*3.25)))
        self.finish()

    def fraunhofer_snr(self):
        self.start_slide('Fraunhofer SNR estimates');self.play(FadeIn(self.title('Fraunhofer tests: conditional thermal budget at 1 m')))
        n=self.h['three_mm/chirps'][:];a=self.h['parallel/ideal_1m_snr_db'][:];b=self.h['perpendicular/ideal_1m_snr_db'][:]
        self.number_table(['Chirps','6606: 5 mm, 78 GHz','6607: 6 mm, 79.5 GHz'],[[k,f'{aa:.1f} dB',f'{bb:.1f} dB'] for k,aa,bb in zip(n,a,b)],.25)
        self.play(FadeIn(self.prose('Maximum power; 6 dBi each antenna; 10 000 K; exact sphere RCS',25).move_to(DOWN*2.2)))
        self.play(FadeIn(self.prose('Budget estimates: geometry, pointing, losses and clutter affect actual SNR',24,MUTED).move_to(DOWN*2.9)))
        self.finish()

    def measured(self):
        self.start_slide('Measured undecoded voltage');self.play(FadeIn(self.title('Measured data: stationary echoes, then disturbed clutter')))
        image=ImageMobject(str(CAPTURE/'output/raw_fast_time_6607_frames_34_35_36_37.png')).scale_to_fit_height(5.3).move_to(UP*.1)
        self.play(FadeIn(image))
        self.play(FadeIn(self.prose('Undecoded real component of complex I/Q, in ADC counts',25).move_to(DOWN*2.85)))
        self.play(FadeIn(self.prose('Early steady structure: vessel walls and stationary chamber objects',24,MUTED).move_to(DOWN*3.35)))
        self.finish()

    def subtracted(self):
        self.start_slide('Quiet background subtraction');self.play(FadeIn(self.title('Remove the quiet chamber background, before decoding')))
        image=ImageMobject(str(CAPTURE/'output/background_subtracted_real_6607_frames_34_35_36_37.png')).scale_to_fit_height(5.25).move_to(UP*.15)
        self.play(FadeIn(image))
        self.play(Write(self.equation(r'b(u)=\langle z_n(u)\rangle_{n<1624},\qquad r_n(u)=z_n(u)-b(u)',size=32).move_to(DOWN*2.9)))
        self.play(FadeIn(self.prose('Same complex background subtracted from quiet and disturbed data',23,MUTED).move_to(DOWN*3.4)))
        self.finish()

    def chirp_spectra(self,shot,frames,subtract):
        suffix=f'{shot}_frames_'+'_'.join(map(str,frames))
        filename=('background_subtracted_spectra_' if subtract else 'chirp_spectra_')+suffix+'.png'
        name=f'Chirp spectra {shot} '+('background subtracted' if subtract else 'raw')
        script='plot_background_subtracted_views.py' if subtract else 'plot_chirp_spectra.py'
        base.NOTES[name]=f'Shot {shot}, RX0, frames {frames}: one complex Hann-windowed FFT128 per chirp, fs12.5MS/s, 97.65625kHz frequency bins. Signed stored-I/Q beat frequency, not range or carrier Doppler alone. Secondary x-axis uses the assumed13.71us chirp period and omits frame gaps. Both before/after plots use -70..0dB relative to the same raw four-frame maximum within this shot. '+('Fixed complex quiet mean removed BEFORE FFT from all chirps, including disturbed data; subtraction is in voltage, not spectral power. ' if subtract else 'No background subtraction. ')+f'Source: capture_analysis/{script}; output/{filename}; matching HDF5 products preserve provenance.'
        self.start_slide(name)
        self.play(FadeIn(self.title(f'Shot {shot}: '+('spectra after quiet-background subtraction' if subtract else 'spectrum of every recorded chirp'))))
        image=ImageMobject(str(CAPTURE/'output'/filename)).scale_to_fit_height(5.25).move_to(UP*.1)
        self.play(FadeIn(image))
        self.play(FadeIn(self.prose('One spectrum per complex chirp; signed beat frequency, not range or Doppler alone',24).move_to(DOWN*2.85)))
        self.play(FadeIn(self.prose('Same colour reference as the raw spectra: compare before and after subtraction' if subtract else 'Hann window; 97.7 kHz bins; time in seconds above; frame gaps omitted',23,MUTED).move_to(DOWN*3.4)))
        self.finish()

    def clutter(self):
        self.start_slide('Disturbance and self noise');self.play(FadeIn(self.title('After the shot: changing chamber echoes can mask the pellet')))
        image=ImageMobject(str(CAPTURE/'output/quiet_subtracted_6607_transition.png')).scale_to_fit_height(4.8).move_to(UP*.45);self.play(FadeIn(image))
        for i,line in enumerate(['Chamber motion changes strong wall echoes: structured “self noise”.',
                                 'A fixed mean removes stationary clutter; changing clutter remains.',
                                 'Thermal SNR alone does not establish that a pellet echo is detectable.']):
            self.play(FadeIn(self.prose(line,25,FG if i<2 else GOLD).move_to(DOWN*(2.3+.55*i))))
        self.finish()

    def simulator(self):
        self.text_slide('Interactive simulator','Explore the FMCW beat signal',[
            'Adjust the target trajectory and play one to ten consecutive chirps.',
            'Follow the target position and the active chirp, including sample gaps.',
            'Compare complex I/Q, instantaneous frequency and the spectrogram.',
            'https://juha.no/fmcw/sim/',
            'Open the simulator using the clickable link at the bottom left.'])

    def matched_noise_definition(self):
        name='Matched energy and noise reference'
        base.NOTES[name]='The previous plot divided matched output by total event residual energy and showed a percentage of explained energy. New rho=|q.H z|^2 / mean_k |q.H n_k|^2, unit-norm Eq25 templates. Denominator estimated separately for every template and I/Q orientation from80disjoint four-chirp quiet reference blocks, accounting for noise colour and within-block correlations. Background complex mean estimated using frame34chirps0–255 only; training256–895, held-out896–1535, event1608–1623 disjoint. Observed output includes signal and noise: under stationary proper complex Gaussian noise, fixed-filter linear mean rho=1 (0dB); with a matched signal, expectation=1+outputSNR. Mean of log(rho) is not0dB. Grid search/maximization raises peak noise response. No universal13dB noise peak or detection threshold. ADC counts cannot by themselves supply calibrated thermal noise at antenna input. Normalization uses total quiet background, including residual clutter, not just thermal noise. This is noise-normalized phase matching, not optimal covariance whitening. Source: capture_analysis/plot_noise_normalized_four_chirps.py; https://www.mathworks.com/help/phased/ug/detection-of-known-signals-with-coherent-and-non-coherent-receivers.html'
        self.start_slide(name)
        self.play(FadeIn(self.title('Matched output energy divided by noise energy')))
        self.play(Write(self.equation(r'\rho(\theta)=\frac{|q_\theta^H z|^2}{\left\langle |q_\theta^H n_k|^2\right\rangle_k},\qquad \rho_{\rm dB}=10\log_{10}\rho',size=43,color=GOLD).move_to(UP*1.95)))
        lines=[('Numerator: energy captured by a trial trajectory filter.',FG),
               ('Denominator: mean output energy for that filter on separate quiet trains.',FG),
               ('One fixed filter: mean noise response is 1, or 0 dB.',RECEIVE),
               ('Searching many filters raises the largest noise peak; 13 dB is not universal.',GOLD),
               ('Here the noise reference includes residual chamber clutter.',MUTED)]
        for i,(line,color) in enumerate(lines):
            self.play(FadeIn(self.prose(line,27,color).move_to(UP*(.55-.75*i))),run_time=.5)
        self.finish()

    def four_chirp_search(self):
        path=CAPTURE/'output/noise_normalized_four_chirps_6607.h5'
        with h5py.File(path) as h:
            summary=json.loads(h.attrs['summary_json'])
        for report in summary['event_reports']:
            start=report['start_chirp'];stop=report['end_chirp']
            name=f'Four chirp noise-normalized match {start}'
            base.NOTES[name]='Shot6607 RX0 frame34, eight fixed neighbouring four-chirp trains1608–1639: four before and four after data-selected disturbance boundary1624. Upper context is the undecoded real component after the SAME complex quiet row-mean cancellation as the filter, in ADC counts. Quiet mean fromchirps0–255 is subtracted from every context chirp, including the disturbed interval. Context extends1560–1999, about5.15ms into the disturbed interval. Gold boundary lines and tint mark analyzed four-chirp train; dashed line is data-selected disturbance boundary, not an independently recorded trigger. Filtering uses complex I/Q minus independent early baseline0–255. Noise reference:80four-chirp blocks every8chirps in256–895. Unit-norm Eq25 templates, per-template/per-IQ noise denominator; colour10log10 of observed matched energy / mean quiet-background energy, common colour limits across all eight trains. Max over y0 for(v0,x0), max over x0 for(v0,y0), IQ orientation also maximized. x0[-1.2,1.2]m step.05m, v0[-7,7]km/s step.1km/s, y0[.05,1]m step.05m. x0 is along-track position, not slant range, relative to t=0 of each train. Exact(x0,v0)mirror ambiguity. T_int40.96us, elapsed sample span51.37us. Period13.71us inferred from reader. Receiver approximation HPFs175/350kHz, absIF<10MHz, at least80percentvalidsamples. These are template matches, not confirmed projectile parameters or calibrated thermal SNR. '+json.dumps(report)+' Script: capture_analysis/plot_noise_normalized_four_chirps.py; HDF5: '+str(path)
            self.start_slide(name)
            self.play(FadeIn(self.title(f'Four neighbouring chirps: {start}--{stop}')))
            figure=ImageMobject(str(CAPTURE/f'output/noise_match_6607_chirps_{start}_{stop}.png')).scale_to_fit_height(5.65).move_to(UP*.05)
            self.play(FadeIn(figure))
            self.play(FadeIn(self.prose('Colour: matched energy / quiet-background noise; maximum over the omitted coordinate',24,GOLD).move_to(DOWN*3.2)))
            self.finish()

    def matched_noise_controls(self):
        path=CAPTURE/'output/noise_normalized_four_chirps_6607.h5'
        with h5py.File(path) as h:
            summary=json.loads(h.attrs['summary_json'])
        name='Quiet-data matched filter controls'
        base.NOTES[name]='Same entire finite Eq25 bank applied to80held-out real quiet four-chirp trains inframe34chirps896–1535. Each train maximized over x0,v0,y0,andstoredIQorientation. Separate noise reference256–895 and baseline0–255.128proper-complex Gaussian surrogate trains use uncentred second moment of the measured training residual, including its spectral colour and cross-chirp correlation. Gaussianity and stationarity not asserted. EmpiricalCDF of searched peaks, not a calibrated rare-false-alarm threshold;80realblocks and128simulationdraws do not support very small Pfa claims. Held-out median10.19dB,max13.56dB; Pre-disturbance trains10.73,12.19,11.79,12.71dB are shown in this control comparison; the four disturbed trains are not included in its vertical reference lines. Pre-disturbance peaks overlap the quiet background. No confirmed pellet echo. Source script capture_analysis/plot_noise_normalized_four_chirps.py; '+json.dumps(summary)
        self.start_slide(name)
        self.play(FadeIn(self.title('How large are peaks from quiet background?')))
        figure=ImageMobject(str(CAPTURE/'output/four_chirp_noise_controls_6607.png')).scale_to_fit_height(5.45).move_to(UP*.1)
        self.play(FadeIn(figure))
        self.play(FadeIn(self.prose('Separate quiet trains: median peak 10.2 dB; largest observed peak 13.6 dB',26).move_to(DOWN*2.85)))
        self.play(FadeIn(self.prose('The pre-disturbance matches overlap the quiet-background peaks',25,GOLD).move_to(DOWN*3.4)))
        self.finish()
