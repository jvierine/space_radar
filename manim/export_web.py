"""Export browser decks with content-addressed media and endpoint validation."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
import tempfile

import cv2
import numpy as np
from PIL import Image, ImageDraw

DECKS = {
    "RadialFFTSearch": ("radial-fft", "Radial FFT matched filtering and laboratory GUI"),
    "FMCWSpaceRadar": ("fmcw-space-radar", "FMCW space radar: power, scattering, integration and chamber data"),
    "RadarEquationNoise": ("radar-equation-noise", "Radar equation, thermal noise and coherent integration"),
}


def credit_author(page):
    """Visible deck credit, outside the animation and navigation controls."""
    page = page.replace("</head>", '<meta name="author" content="Juha Vierinen">\n</head>')
    credit = ('<span id="deck-author" style="position:fixed;right:110px;bottom:14px;'
              'z-index:1000;background:white;color:#176BB0;padding:6px 10px;'
              'border-radius:4px;font:14px Arial">Juha Vierinen</span>')
    return page.replace("</body>", credit + "\n</body>")


def export(scene, root):
    slug, title = DECKS[scene]
    dest = root / slug
    dest.mkdir(parents=True, exist_ok=True)
    slide_config = json.loads(Path(f"slides/{scene}.json").read_text())
    render_folders = {Path(slide["file"]).resolve().parent for slide in slide_config["slides"]}
    if len(render_folders) != 1:
        raise RuntimeError(
            f"{scene}: mixed scene renders are unsafe. A following segment can "
            "still contain a fade-out of obsolete content. Re-render the entire "
            "deck in execution order before exporting."
        )
    frames = []
    for number, slide in enumerate(json.loads(Path(f"slides/{scene}.json").read_text())["slides"], 1):
        cap = cv2.VideoCapture(slide["file"])
        cap.set(cv2.CAP_PROP_POS_FRAMES, max(0, int(cap.get(cv2.CAP_PROP_FRAME_COUNT)) - 1))
        ok, frame = cap.read()
        cap.release()
        if not ok:
            raise RuntimeError(f"Unreadable slide: {scene} {number}")
        fraction = (np.max(np.abs(frame.astype(np.int16) - frame[0, 0].astype(np.int16)), axis=2) > 18).mean()
        if fraction < 0.01:
            raise RuntimeError(f"Blank slide ending: {scene} {number}")
        picture = Image.fromarray(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB))
        picture.save(dest / f"qa-{number:02d}.png")
        picture.thumbnail((640, 360))
        frames.append(picture)
    contact = Image.new("RGB", (1280, ((len(frames) + 1) // 2) * 388), "#253040")
    draw = ImageDraw.Draw(contact)
    for i, frame in enumerate(frames):
        x, y = (i % 2) * 640, (i // 2) * 388
        contact.paste(frame, (x, y))
        draw.text((x + 10, y + 365), f"Slide {i + 1}", fill="white")
    contact.save(dest / "qa-contact.jpg")
    # Isolated scene renders can reuse Manim's basename across different folders.
    # The HTML converter flattens those folders, so disambiguate BEFORE copying.
    with tempfile.TemporaryDirectory(prefix="manim-web-export-") as staging:
        staging = Path(staging)
        config = json.loads(Path(f"slides/{scene}.json").read_text())
        for slide in config["slides"]:
            for key in ("file", "rev_file"):
                source = Path(slide[key])
                digest = hashlib.sha256(source.read_bytes()).hexdigest()
                target = staging / (digest + ".mp4")
                if not target.exists():
                    shutil.copyfile(source, target)
                slide[key] = str(target)
        (staging / f"{scene}.json").write_text(json.dumps(config))
        subprocess.run(["manim-slides", "convert", "--folder", str(staging), "--offline", scene, str(dest / "index.html")], check=True)
    page = (dest / "index.html").read_text().replace("<title>Manim Slides</title>", f"<title>{title}</title>")
    for media in (dest / "index_assets").glob("*.mp4"):
        new_name = hashlib.sha256(media.read_bytes()).hexdigest()[:20] + ".mp4"
        page = page.replace("index_assets/" + media.name, "index_assets/" + new_name)
        media.rename(media.with_name(new_name))
    guard = Path("final_slide_guard.js").read_bytes()
    guard_name = "navigation-" + hashlib.sha256(guard).hexdigest()[:12] + ".js"
    (dest / "index_assets" / guard_name).write_bytes(guard)
    page = page.replace("</body>", f'<script src="index_assets/{guard_name}"></script>\n</body>')
    if scene == "RadialFFTSearch":
        page = page.replace("</body>", '<a href="https://juha.no/fmcw/lab/" target="_blank" rel="noopener" style="position:fixed;left:18px;bottom:14px;z-index:1000;background:white;color:#176BB0;padding:8px 12px;border:1px solid #176BB0;border-radius:6px;font:16px Arial;text-decoration:none">Open laboratory ↗</a>\n</body>')
    (dest / "index.html").write_text(credit_author(page))
    print(f"{scene}: {len(frames)} nonblank slide endings; exported to {dest}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("scenes", nargs="+", choices=list(DECKS))
    args = parser.parse_args()
    for scene in args.scenes:
        export(scene, args.output)
