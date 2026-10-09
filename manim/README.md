# FMCW space-radar lecture

Run the commands in the root README from this directory. `fmcw_space_radar.py`
contains the complete 47-slide deck; `radar_equation_noise.py` supplies shared
styles, slide boundaries and a reusable introductory scene. The unchanged
TeX-backed `Text` class was extracted from `jvierine/manimations` into
`planck_to_ktb.py`, so the lecture does not depend on unrelated teaching decks.

`space_radar_assets.py` regenerates `assets/fmcw_lecture.h5` from the bundled
synthetic `../capture_analysis/search_experiment.h5`. It checks exact complex
waveform agreement with the reference Python model before writing the SNR,
conducting-sphere scattering and three-chirp matched-filter examples.

## Interpretation

- Maximum TX power: 12 dBm per TX; one TX and one RX. Antenna gains 6–10 dBi
  toward the target and system temperature 10000 K remain budget assumptions.
- x0 is the radar's longitudinal offset relative to target x=0 at the start
  of the selected train. y0 is perpendicular distance; initial slant range
  is sqrt(x0^2+y0^2). Signed v0 is along-track velocity, not radial velocity.
- 128 complex samples at 12.5 MS/s provide 10.24 microseconds per chirp.
  The 13.71 microsecond repetition interval is inferred from reader timing;
  elapsed time plots omit frame gaps. T_int denotes coherently used sample
  duration, excluding unsampled gaps. The constant-amplitude, two-sided
  complex-baseband noise-equivalent bandwidth is 1/T_int.
- The FDTD GIF is from `sabmod`: `figures/maarsy_rcs/core_dense_dx32_e_field_zoom.gif`,
  generated with `animate_e_field.py`, `maarsy_rcs.py` and `fdtd/src/main.rs`.
  It illustrates meteor-plasma total/scattered Ez at 53.5 MHz. Its unresolved
  dense core is a visualization, not a calibrated RCS calculation or an
  IWR1843 pellet simulation. The source remains in the sabmod project.
- Fraunhofer chamber geometry is schematic, not surveyed. Radar 1 is at the
  left end looking toward the impact target. The one-metre useful trajectory
  window is an assumption, not a measured beam footprint.
- The measured four-chirp searches use the same complex quiet mean from
  chirps 0–255 for both the upper context and the filter. Eight consecutive
  trains cover 1608–1639 around the data-selected disturbance boundary 1624.
  Upper contexts extend through chirp 1999 without FFT/decoding.
- Search colours are observed matched energy divided by independent mean
  quiet-background matched energy, in dB. The denominator includes residual
  clutter. This is not calibrated thermal SNR or a confirmed pellet detection.
  Both 2D surfaces maximize over the omitted coordinate and I/Q orientation.
  The mirrored pair (x0,v0) and (-x0,-v0) is indistinguishable in this geometry.

Presenter notes contain further assumptions and primary-source links. The
figure manifest in `../output/figure_provenance.json` maps each bundled plot
to its generating script. Rendered movies and HTML exports are generated
locally and excluded from Git; published browser assets remain at juha.no.

## Radial FFT search and GUI

`radial_fft_search.py` is a separate 12-slide white-background deck. It explains
the three motion parameters, straight-line geometry, phase alignment,
correction-plus-FFT search, per-receiver noise estimation, explained matched
powers and map projections, and graphical selections/results/export. All illustrated
trajectories and phasors are synthetic. The slides link to the laboratory; technical derivations remain in the memos.

```bash
conda run -n base python radial_fft_assets.py
conda run --no-capture-output -n base manim-slides render --disable_caching -r 1920,1080 --fps 30 radial_fft_search.py RadialFFTSearch
conda run -n base python export_web.py qa_fmcw_web RadialFFTSearch
```

From the repository root, `conda run -n base python tools/radial_complexity.py`
regenerates `manim/assets/radial_complexity.h5` for memo 10. The shipped planner
counts are independently checked by `node tools/validate_complexity_counts.mjs`.
