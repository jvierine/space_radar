Use the base conda environment for Python commands (`conda run -n base python ...`).
Use HDF5, not CSV, for scientific data products.

Keep script provenance visible by default in scientific figures and slides;
provide a toggle when embedding scientific figures in an article or memo.

For Manim slides, pause on completed content before fading it out. Preserve
the deferred cleanup in start_slide/clear_slide and final_slide_guard.js.
Rebuild the full affected scene rather than splicing clips into an old manifest.
Inspect every exported slide's final decoded frame. Use the shared TeX-backed
Text renderer in manim/planck_to_ktb.py, imported after `from manim import *`.

The FMCW web destinations are https://juha.no/fmcw/sim/ and
https://juha.no/fmcw/slides/. Keep them unlisted from the /space/ course index.
Verify natural browser playback and final-slide retention after deployment.
Do not publish raw recordings, intermediate capture arrays or generated videos
as source-code assets without an explicit request.
