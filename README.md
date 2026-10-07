# FMCW space radar

Source for the [interactive simulator](https://juha.no/fmcw/sim/) and
[Manim lecture](https://juha.no/fmcw/slides/), plus their Python signal model
and scientific figure scripts.

## Interactive simulator

From the repository root:

```sh
conda run --no-capture-output -n base python -m http.server 18765 --bind 127.0.0.1 --directory web
```

Open **http://127.0.0.1:18765/sim/**. No npm installation or build is needed;
Three.js and its MIT license are bundled in `web/sim/vendor/`.

The simulator shows target position, 1–10 chirps with ADC/tail/idle timing,
analytical instantaneous frequency, complex beat voltage after analogue
chirp downconversion, and a spectrogram. Playback synchronizes the target,
active chirp and plot cursors. It uses the memo's Equation 25 with retarded
target geometry; the voltage scale is relative, not calibrated ADC volts.
See [model details](web/sim/README.md).

Validate the JavaScript model against Python (requires Node.js):

```sh
conda run -n base python web/sim/validate_model.py
conda run -n base python -m unittest discover -s capture_analysis -p 'test_*.py'
```

## Manim slides

The complete **47-slide** scene is
[`manim/fmcw_space_radar.py`](manim/fmcw_space_radar.py), class
`FMCWSpaceRadar`. The bundled lecture inputs allow rendering without the
original capture disk or another repository. Use the base conda environment.
Python dependencies are recorded in `requirements.txt`; Manim additionally
needs FFmpeg, LaTeX and `dvisvgm` on the system path.

```sh
cd manim
conda run -n base python space_radar_assets.py
conda run --no-capture-output -n base manim-slides render --disable_caching -r 1920,1080 --fps 30 fmcw_space_radar.py FMCWSpaceRadar
conda run -n base python export_web.py qa_fmcw_web FMCWSpaceRadar
conda run -n base python build_radar_calculator.py
conda run --no-capture-output -n base python -m http.server 18766 --bind 127.0.0.1 --directory qa_fmcw_web
```

Open **http://127.0.0.1:18766/fmcw-space-radar/**. The exporter checks the
last decoded frame of every slide, saves a contact sheet, hashes video assets,
and installs the final-slide navigation guard. Source footers are shown by
default; set `SHOW_PROVENANCE=0` to hide them.

The deck derives the radar equation, scattering and thermal noise; explains
coherent integration; demonstrates the simulated beat voltage and matched
filters; and shows measured voltage, per-chirp spectra and eight neighbouring
four-chirp template searches across the disturbance boundary. The 2D search
plots take the **maximum** over the omitted coordinate and I/Q orientation.
See [lecture inputs and interpretation](manim/README.md).

## Data and provenance

- `capture_analysis/` contains the Equation 25 reference model, scattering
  solver, simulation and measurement-analysis scripts, and a synthetic HDF5
  experiment used to regenerate the lecture's simulated signals.
- `output/` contains only the figures already used in the lecture, their
  SHA-256/generator manifest, and a small HDF5 summary for rendering. The
  summary is explicitly marked as presentation metadata; it does **not**
  contain the full matched-filter search cube or measured I/Q samples.
- `manim/assets/` contains the generated lecture HDF5 and the existing FDTD
  demonstration GIF. The GIF depicts meteor plasma at 53.5 MHz, not the metal
  pellet at 78 GHz; see its provenance in the lecture notes.
- Original Fraunhofer recordings and intermediate capture products are not
  included. For reanalysis, set `FMCW_CAPTURE_DIR` to the directory containing
  the original recordings. Some analysis scripts also require earlier HDF5
  products in `output/`; see [analysis instructions](capture_analysis/README.md).

Existing `read_subframe_res.py` is retained as the original exploratory reader;
it still expects its original dataset path. The portable simulator and lecture
use the sources documented above.
