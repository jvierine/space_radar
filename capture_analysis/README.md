# Scientific reference model and analysis

Analysis programs by **Juha Vierinen**.

Use `conda run -n base python ...` from the repository root. Data products are
HDF5. `analyze_captures.py` defines the exact retarded-range Equation 25 echo,
its instantaneous frequency, receiver approximation and PEC-sphere solver.
The bundled `search_experiment.h5` is synthetic, including complex simulated
signals and noise, not original measured recordings.

Independent numerical checks:

```sh
conda run -n base python -m unittest discover -s capture_analysis -p 'test_*.py'
conda run -n base python web/sim/validate_model.py
```

The full lecture can be rendered without re-running capture analysis. To
reanalyse measured data, provide original files via `FMCW_CAPTURE_DIR`
(default `data/raw/`). The analysis scripts retain the established shot/file
names and HDF5 schemas. No raw measurements or intermediate capture products
are distributed in this source repository.

Typical analysis chain, with the original inputs present:

1. `analyze_captures.py`: profiles and feasibility products.
2. `search_simulated_echoes.py`: synthetic echo/template experiments.
3. `plot_raw_fast_time.py`, `plot_chirp_spectra.py`,
   `plot_background_subtracted_views.py`: voltage/spectrum products. Run each
   script with `--help` where supported and inspect its input paths.
4. `plot_raw_receiver_iq.py`, `plot_quiet_subtracted.py`: receiver-profile and
   quiet-subtracted products needed by subsequent diagnostics.
5. `plot_noise_normalized_four_chirps.py`: reads the per-chirp spectrum HDF5
   containing full complex counts, `raw_receiver_iq.h5`, and the background
   subtraction HDF5; generates eight four-chirp searches and control plots.

Historical wider searches and diagnostics are included as source, with their
existing input/output relationships. Signed-velocity GPU/CPU diagnostics
optionally require PyTorch. Earlier analysis products are not automatically
recreated by rendering the deck.

`../output/noise_normalized_four_chirps_6607.h5` contains only the summary used
to render slides, marked `presentation_metadata_only`. It is not an analysis
input or a substitute for the full search product. Re-running the search
with the original inputs writes the complete scientific HDF5 in its place.
