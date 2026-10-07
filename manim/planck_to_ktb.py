"""Shared TeX-backed Text renderer extracted unchanged from jvierine/manimations.
Only the renderer needed by this radar deck is vendored here.
"""
from manim import *


class Text(Tex):
    """Typeset prose with LaTeX; bypass Pango/SVG glyph-position corruption.

    Keep the deck's serif appearance and TeX's word spacing/kerning. Escape
    ordinary prose before typesetting; mathematical content uses MathTex.
    """

    def __init__(self, text, font_size=48, color=WHITE, weight=None, t2c=None, **kwargs):
        escapes = {
            "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
            "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
            "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
            "–": "--", "—": "---", "\n": r"\\",
        }
        def escape(value):
            return "".join(escapes.get(char, char) for char in value)
        prose = escape(text)
        for phrase, shade in (t2c or {}).items():
            rgb = ManimColor(shade).to_hex().lstrip("#")
            prose = prose.replace(escape(phrase), r"\textcolor[HTML]{" + rgb + "}{" + escape(phrase) + "}")
        if weight in (BOLD, SEMIBOLD):
            prose = r"\textbf{" + prose + "}"
        template = TexTemplate()
        template.add_to_preamble(r"\usepackage{xcolor}")
        kwargs.setdefault("tex_template", template)
        super().__init__(prose, font_size=font_size, color=color, **kwargs)
