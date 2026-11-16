#!/usr/bin/env python3
"""Common rendering layer for plant-virus-screen graphic summaries.

One Figure primitive model renders the same drawing to:

- a self-contained HTML page with an inline SVG (for browsing), and
- a native vector PDF written directly by this module (text, lines and
  shapes stay sharp and scalable at any zoom; no LibreOffice/soffice).

The nt and the protein visualizers build their figures from these shared
primitives and keep only their layout-specific code, so the duplicated
SVG-string/PDF-export logic from the earlier per-script implementations
lives here exactly once.

Python 3.6+ compatible, standard library only.

write_output_bundle keeps the historical signature
``(out_path, html_text, write_svg=True, write_pdf=True)`` for backward
compatibility; standalone SVG output has been replaced by native vector
PDF, so ``write_svg`` is accepted but ignored.
"""

import colorsys
import html
import io
import math
import os
from pathlib import Path

# CSS px (96 dpi) -> PDF points (72 dpi).
PT_PER_PX = 0.75

# ---------------------------------------------------------------------------
# Helvetica metrics (Adobe core14 AFM advance widths, chars 32..126,
# units per 1000). Used for text anchoring in the PDF backend; Arial used
# by the SVG backend has nearly identical metrics.
# ---------------------------------------------------------------------------

_HELVETICA = (
    278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278,
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556, 1015,
    667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778, 667, 778,
    722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556, 333,
    556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556, 556,
    556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584,
)

_HELVETICA_BOLD = (
    278, 333, 474, 556, 556, 889, 722, 238, 333, 333, 389, 584, 278, 333, 278, 278,
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 333, 333, 584, 584, 584, 611, 975,
    722, 722, 722, 722, 667, 611, 778, 722, 278, 556, 722, 611, 833, 722, 778, 667, 778,
    722, 667, 611, 722, 667, 944, 667, 667, 611, 333, 278, 333, 584, 556, 333,
    556, 611, 556, 611, 556, 333, 611, 611, 278, 278, 556, 278, 889, 611, 611, 611,
    611, 389, 556, 333, 611, 556, 778, 556, 556, 500, 389, 280, 389, 584,
)

_FALLBACK_WIDTH = 0.55

# ---------------------------------------------------------------------------
# CJK font support for the PDF backend
#
# Text containing non-latin-1 characters is rendered through an embedded
# subset of an open-licence CJK TrueType font (default search order below,
# override with the PVS_CJK_FONT environment variable). Embedding the
# subset guarantees correct display in every PDF viewer; when no usable
# font + fontTools installation is found we fall back to a non-embedded
# Type0 reference (STSong-Light) that viewers substitute themselves.
# ---------------------------------------------------------------------------

_CJK_FONT_CANDIDATES = [
    "/usr/share/fonts/wqy-microhei/wqy-microhei.ttc",
    "/usr/share/fonts/wqy-zenhei/wqy-zenhei.ttc",
    "/usr/share/fonts/cjkuni-uming/uming.ttc",
]

_CJK_METRICS_CACHE = {}


def _find_cjk_font_path():
    override = os.environ.get("PVS_CJK_FONT")
    if override and Path(override).is_file():
        return override
    for candidate in _CJK_FONT_CANDIDATES:
        if Path(candidate).is_file():
            return candidate
    return None


def _cjk_metrics():
    """Unicode -> (glyph id, advance width, unitsPerEm) map, cached.

    Returns None when fontTools or a CJK font is unavailable.
    """
    if _CJK_METRICS_CACHE:
        return _CJK_METRICS_CACHE.get("metrics")
    metrics = None
    path = _find_cjk_font_path()
    if path:
        try:
            from fontTools.ttLib import TTFont

            # fontNumber=0 selects the first font of a TTC collection and is
            # ignored for plain TTF/OTF files.
            font = TTFont(path, fontNumber=0)
            cmap = font.getBestCmap()
            hmtx = font["hmtx"]
            upem = font["head"].unitsPerEm
            gid = {}
            adv = {}
            for code, glyph_name in cmap.items():
                gid[code] = font.getGlyphID(glyph_name)
                adv[code] = hmtx[glyph_name][0]
            metrics = {"path": path, "gid": gid, "adv": adv, "upem": upem}
        except Exception:
            metrics = None
    _CJK_METRICS_CACHE["metrics"] = metrics
    return metrics


def _subset_cjk_font(used_chars):
    """Subset the CJK font to ``used_chars``; returns (ttf_bytes, metrics).

    A fresh TTFont is opened every call so the cached metrics font object
    is never mutated; retain_gids keeps glyph numbering stable, which is
    what the embedded Identity-H text codes rely on.
    """
    metrics = _cjk_metrics()
    from fontTools import subset as _subset
    from fontTools.ttLib import TTFont

    font = TTFont(metrics["path"], fontNumber=0)
    unicodes = sorted(ord(ch) for ch in used_chars if ord(ch) in metrics["gid"])
    options = _subset.Options()
    options.retain_gids = True
    subsetter = _subset.Subsetter(options)
    subsetter.populate(unicodes=unicodes)
    subsetter.subset(font)
    buffer = io.BytesIO()
    font.save(buffer)
    return buffer.getvalue(), metrics


def _has_nonlatin(s):
    try:
        s.encode("latin-1")
    except UnicodeEncodeError:
        return True
    return False


def _char_width(ch, bold):
    code = ord(ch)
    if 32 <= code <= 126:
        table = _HELVETICA_BOLD if bold else _HELVETICA
        return table[code - 32] / 1000.0
    return _FALLBACK_WIDTH


def text_width(s, size, bold=False):
    """Approximate rendered width (in figure px) of ``s`` at font ``size``.

    Latin-1 text uses Helvetica metrics; text containing CJK glyphs uses
    the loaded CJK font advances when available, otherwise a
    full-width/half-width approximation.
    """
    if not s:
        return 0.0
    try:
        s.encode("latin-1")
    except UnicodeEncodeError:
        metrics = _cjk_metrics()
        if metrics:
            total = 0.0
            for ch in s:
                advance = metrics["adv"].get(ord(ch))
                if advance is None:
                    total += 0.5 if ord(ch) < 0x2E80 else 1.0
                else:
                    total += advance / float(metrics["upem"])
            return total * float(size)
        width = 0.0
        for ch in s:
            width += 0.5 if ord(ch) < 0x2E80 else 1.0
        return width * float(size)
    return sum(_char_width(ch, bold) for ch in s) * float(size)


# ---------------------------------------------------------------------------
# Color & text layout helpers (Nature-figure sequential palette & wrapping)
# ---------------------------------------------------------------------------

_NATURE_IDENTITY_STOPS = (
    (0.00, (156, 202, 214)),  # #9CCAD6 soft ice-teal (high luminance)
    (0.35, (66, 148, 158)),   # #42949E Nature muted teal
    (0.70, (55, 117, 186)),   # #3775BA Nature secondary blue
    (1.00, (15, 77, 146)),    # #0F4D92 Nature primary deep navy (low luminance)
)


def hex_to_rgb01(color):
    color = (color or "#000000").strip()
    if color.startswith("#"):
        color = color[1:]
    if len(color) == 3:
        color = "".join(ch * 2 for ch in color)
    value = int(color, 16)
    return ((value >> 16) & 255) / 255.0, ((value >> 8) & 255) / 255.0, (value & 255) / 255.0


def pident_to_color(pident, low, high):
    """Map an identity percentage onto a colorblind-safe Nature sequential ramp.

    Uses a monotonically darkening ice-teal -> teal -> blue -> deep-navy ramp
    so both colour and grayscale print preserve the exact identity ordering.
    """
    if high <= low:
        ratio = 0.5
    else:
        ratio = (pident - low) / (high - low)
    ratio = max(0.0, min(1.0, float(ratio)))
    for idx in range(len(_NATURE_IDENTITY_STOPS) - 1):
        t0, rgb0 = _NATURE_IDENTITY_STOPS[idx]
        t1, rgb1 = _NATURE_IDENTITY_STOPS[idx + 1]
        if ratio <= t1 or idx == len(_NATURE_IDENTITY_STOPS) - 2:
            span = max(1e-9, t1 - t0)
            local = max(0.0, min(1.0, (ratio - t0) / span))
            r = int(round(rgb0[0] + (rgb1[0] - rgb0[0]) * local))
            g = int(round(rgb0[1] + (rgb1[1] - rgb0[1]) * local))
            b = int(round(rgb0[2] + (rgb1[2] - rgb0[2]) * local))
            return "#{:02X}{:02X}{:02X}".format(r, g, b)
    r, g, b = _NATURE_IDENTITY_STOPS[-1][1]
    return "#{:02X}{:02X}{:02X}".format(r, g, b)


def wrap_text(text, max_width, size=11, bold=False):
    """Wrap text to fit within ``max_width`` figure px using Helvetica/CJK metrics."""
    lines = []
    remaining = str(text or "").strip()
    max_width = max(24.0, float(max_width))
    while remaining:
        if text_width(remaining, size, bold=bold) <= max_width:
            lines.append(remaining)
            break
        cut = 1
        while cut < len(remaining) and text_width(remaining[: cut + 1], size, bold=bold) <= max_width:
            cut += 1
        break_pos = -1
        for sep in (" ", "|", "_", "-", "/", ","):
            pos = remaining.rfind(sep, 0, cut + 1)
            if pos > cut // 2:
                break_pos = pos + (0 if sep == " " else 1)
                break
        if break_pos > 0:
            cut = break_pos
        lines.append(remaining[:cut].rstrip())
        remaining = remaining[cut:].lstrip()
    return lines or [""]


def assign_text_lanes(entries, pad_px=10.0):
    """Greedy interval-packing of horizontal spans ``[span_left, span_right]`` into lanes."""
    lane_last_end = []
    for entry in entries:
        placed = False
        for lane_index, last_end in enumerate(lane_last_end):
            if entry["span_left"] > last_end + pad_px:
                entry["lane"] = lane_index
                lane_last_end[lane_index] = entry["span_right"]
                placed = True
                break
        if not placed:
            entry["lane"] = len(lane_last_end)
            lane_last_end.append(entry["span_right"])
    return len(lane_last_end)


def observed_identity_range(grouped_hits):
    """Identity range across all hits, padded when a single value occurs."""
    values = [hit.pident for hits in grouped_hits.values() for hit in hits]
    if not values:
        return 0.0, 100.0
    low = min(values)
    high = max(values)
    if math.isclose(low, high):
        low = max(0.0, low - 5.0)
        high = min(100.0, high + 5.0)
    return low, high


# ---------------------------------------------------------------------------
# Figure primitive model
# ---------------------------------------------------------------------------


class Figure(object):
    """Ordered list of drawing primitives in figure pixel coordinates.

    y grows downward (SVG convention); the PDF renderer flips it.

    Primitives: text, rect, polyline, gradient_rect. Tooltips are recorded
    as ``title`` and rendered as SVG <title> elements (PDF has no tooltips).
    ``rx`` rounded corners are kept for SVG; the PDF backend draws plain
    rectangles (visually negligible at these radii).
    """

    def __init__(self, width, height):
        self.width = float(width)
        self.height = float(height)
        self.ops = []

    def text(self, x, y, s, size=12, fill="#222222", anchor="start", bold=False, title=None):
        self.ops.append(
            {
                "op": "text",
                "x": float(x),
                "y": float(y),
                "s": str(s),
                "size": float(size),
                "fill": fill,
                "anchor": anchor,
                "bold": bool(bold),
                "title": title,
            }
        )

    def rect(self, x, y, w, h, fill=None, stroke=None, stroke_width=1.0, rx=0.0, title=None):
        self.ops.append(
            {
                "op": "rect",
                "x": float(x),
                "y": float(y),
                "w": float(w),
                "h": float(h),
                "fill": fill,
                "stroke": stroke,
                "stroke_width": float(stroke_width),
                "rx": float(rx),
                "title": title,
            }
        )

    def polyline(self, points, stroke="#111111", stroke_width=1.0, dash=None, title=None):
        self.ops.append(
            {
                "op": "polyline",
                "points": [(float(px), float(py)) for px, py in points],
                "stroke": stroke,
                "stroke_width": float(stroke_width),
                "dash": dash,
                "title": title,
            }
        )

    def polygon(self, points, fill="#bbbbbb", stroke=None, stroke_width=1.0):
        """Closed filled polygon (e.g. area under a curve)."""
        self.ops.append(
            {
                "op": "polygon",
                "points": [(float(px), float(py)) for px, py in points],
                "fill": fill,
                "stroke": stroke,
                "stroke_width": float(stroke_width),
            }
        )

    def circle(self, cx, cy, r, fill="#111111", stroke=None, stroke_width=1.0, title=None):
        """Circle (scatter points); PDF backend approximates with 4 bezier arcs."""
        self.ops.append(
            {
                "op": "circle",
                "cx": float(cx),
                "cy": float(cy),
                "r": float(r),
                "fill": fill,
                "stroke": stroke,
                "stroke_width": float(stroke_width),
                "title": title,
            }
        )

    def gradient_rect(self, x, y, w, h, color_fn, steps=96, stroke="#bbbbbb"):
        """Horizontal gradient bar; color_fn(frac in 0..1) -> '#RRGGBB'."""
        self.ops.append(
            {
                "op": "gradient_rect",
                "x": float(x),
                "y": float(y),
                "w": float(w),
                "h": float(h),
                "color_fn": color_fn,
                "steps": int(steps),
                "stroke": stroke,
            }
        )


class HtmlDoc(str):
    """Full HTML page text that carries its Figure for native PDF export."""

    def __new__(cls, value, figure):
        obj = super(HtmlDoc, cls).__new__(cls, value)
        obj.figure = figure
        return obj


# ---------------------------------------------------------------------------
# SVG / HTML backend
# ---------------------------------------------------------------------------


def _fmt(value):
    return "{:.2f}".format(float(value))


def _title_xml(op):
    if not op.get("title"):
        return ""
    return "<title>{}</title>".format(html.escape(op["title"], quote=True))


def render_svg(fig):
    """Render the figure as a standalone, style-inline <svg> element."""
    gradient_ops = [op for op in fig.ops if op["op"] == "gradient_rect"]
    parts = []
    parts.append(
        "<svg width='{w}' height='{h}' viewBox='0 0 {w} {h}' xmlns='http://www.w3.org/2000/svg' "
        "font-family='Arial,Helvetica,sans-serif' style='background:#ffffff;'>".format(
            w=_fmt(fig.width), h=_fmt(fig.height)
        )
    )
    if gradient_ops:
        parts.append("<defs>")
        for index, op in enumerate(gradient_ops):
            parts.append(
                "<linearGradient id='grad-{i}' x1='0%' y1='0%' x2='100%' y2='0%'>".format(i=index)
            )
            steps = max(2, op["steps"])
            for idx in range(steps):
                frac = float(idx) / max(1, steps - 1)
                parts.append(
                    "<stop offset='{pct:.2f}%' stop-color='{color}'/>".format(
                        pct=frac * 100.0, color=op["color_fn"](frac)
                    )
                )
            parts.append("</linearGradient>")
        parts.append("</defs>")
    for op in fig.ops:
        kind = op["op"]
        if kind == "text":
            attrs = "x='{x}' y='{y}' font-size='{size}' fill='{fill}'".format(
                x=_fmt(op["x"]), y=_fmt(op["y"]), size=_fmt(op["size"]), fill=op["fill"]
            )
            if op["anchor"] != "start":
                attrs += " text-anchor='{}'".format(op["anchor"])
            if op["bold"]:
                attrs += " font-weight='bold'"
            parts.append(
                "<text {}>{}{}</text>".format(attrs, _title_xml(op), html.escape(op["s"]))
            )
        elif kind == "rect":
            attrs = "x='{x}' y='{y}' width='{w}' height='{h}'".format(
                x=_fmt(op["x"]), y=_fmt(op["y"]), w=_fmt(op["w"]), h=_fmt(op["h"])
            )
            if op["rx"]:
                attrs += " rx='{}'".format(_fmt(op["rx"]))
            attrs += " fill='{}'".format(op["fill"] if op["fill"] else "none")
            if op["stroke"]:
                attrs += " stroke='{}' stroke-width='{}'".format(op["stroke"], _fmt(op["stroke_width"]))
            child = _title_xml(op)
            if child:
                parts.append("<rect {}>{}</rect>".format(attrs, child))
            else:
                parts.append("<rect {}/>".format(attrs))
        elif kind == "polyline":
            points = " ".join(
                "{},{}".format(_fmt(px), _fmt(py)) for px, py in op["points"]
            )
            attrs = "points='{points}' fill='none' stroke='{stroke}' stroke-width='{sw}'".format(
                points=points, stroke=op["stroke"], sw=_fmt(op["stroke_width"])
            )
            if op.get("dash"):
                attrs += " stroke-dasharray='{},{}'".format(_fmt(op["dash"][0]), _fmt(op["dash"][1]))
            child = _title_xml(op)
            if child:
                parts.append("<polyline {}>{}</polyline>".format(attrs, child))
            else:
                parts.append("<polyline {}/>".format(attrs))
        elif kind == "polygon":
            points = " ".join(
                "{},{}".format(_fmt(px), _fmt(py)) for px, py in op["points"]
            )
            attrs = "points='{points}' fill='{fill}'".format(points=points, fill=op["fill"])
            if op["stroke"]:
                attrs += " stroke='{}' stroke-width='{}'".format(op["stroke"], _fmt(op["stroke_width"]))
            parts.append("<polygon {}/>".format(attrs))
        elif kind == "circle":
            attrs = "cx='{cx}' cy='{cy}' r='{r}' fill='{fill}'".format(
                cx=_fmt(op["cx"]), cy=_fmt(op["cy"]), r=_fmt(op["r"]), fill=op["fill"]
            )
            if op["stroke"]:
                attrs += " stroke='{}' stroke-width='{}'".format(op["stroke"], _fmt(op["stroke_width"]))
            child = _title_xml(op)
            if child:
                parts.append("<circle {}>{}</circle>".format(attrs, child))
            else:
                parts.append("<circle {}/>".format(attrs))
        elif kind == "gradient_rect":
            index = gradient_ops.index(op)
            attrs = "x='{x}' y='{y}' width='{w}' height='{h}' fill='url(#grad-{i})'".format(
                x=_fmt(op["x"]), y=_fmt(op["y"]), w=_fmt(op["w"]), h=_fmt(op["h"]), i=index
            )
            if op["stroke"]:
                attrs += " stroke='{}'".format(op["stroke"])
            parts.append("<rect {}/>".format(attrs))
    parts.append("</svg>")
    return "\n".join(parts)


def render_html(fig, page_title, h1_text, meta_text):
    """Render a self-contained HTML page wrapping the inline SVG figure."""
    parts = []
    parts.append("<!DOCTYPE html>")
    parts.append("<html lang='en'>")
    parts.append("<head>")
    parts.append("<meta charset='utf-8'>")
    parts.append("<title>{}</title>".format(html.escape(page_title)))
    parts.append(
        "<style>"
        "body{font-family:Arial,Helvetica,sans-serif;margin:24px;color:#222;background:#fff;}"
        "h1{font-size:22px;margin:0 0 8px 0;}"
        ".meta{margin:0 0 18px 0;font-size:14px;color:#444;}"
        "</style>"
    )
    parts.append("</head>")
    parts.append("<body>")
    parts.append("<h1>{}</h1>".format(html.escape(h1_text)))
    parts.append("<div class='meta'>{}</div>".format(meta_text))
    parts.append(render_svg(fig))
    parts.append("</body>")
    parts.append("</html>")
    return HtmlDoc("\n".join(parts), fig)


# ---------------------------------------------------------------------------
# Native vector PDF backend (no external tools)
# ---------------------------------------------------------------------------


def _pdf_text_escape(s):
    s = s.encode("latin-1", "replace").decode("latin-1")
    return s.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _col(color):
    return "{:.4f} {:.4f} {:.4f}".format(*hex_to_rgb01(color))


def render_pdf_stream(fig, cjk_mode="ucs2"):
    """Build the PDF page content stream for the figure (text ops in pt).

    ``cjk_mode``: "none" (no CJK text), "gid" (text codes are glyph ids of
    the embedded subset font, Identity-H), or "ucs2" (fallback: UTF-16BE
    codes for the non-embedded STSong-Light reference).
    """
    k = PT_PER_PX
    height = fig.height
    out = []
    out.append("q")
    # Pure white publication canvas (no outer decorative box border).
    out.append("{} rg".format(_col("#ffffff")))
    out.append("0 0 {w:.2f} {h:.2f} re f".format(w=fig.width * k, h=height * k))
    for op in fig.ops:
        kind = op["op"]
        if kind == "text":
            x = op["x"]
            width = text_width(op["s"], op["size"], op["bold"])
            if op["anchor"] == "middle":
                x -= width / 2.0
            elif op["anchor"] == "end":
                x -= width
            y = height - op["y"]
            cjk = _has_nonlatin(op["s"])
            if cjk:
                font = "/F3"
                if cjk_mode == "gid":
                    metrics = _cjk_metrics()
                    codes = []
                    for ch in op["s"]:
                        glyph = metrics["gid"].get(ord(ch), 0)
                        codes.append("{:04X}".format(glyph))
                    payload = "<{}>".format("".join(codes))
                else:
                    payload = "<{}>".format(op["s"].encode("utf-16-be").hex())
            else:
                font = "/F2" if op["bold"] else "/F1"
                payload = "({})".format(_pdf_text_escape(op["s"]))
            out.append("{} rg".format(_col(op["fill"])))
            out.append("BT")
            out.append("{} {:.2f} Tf".format(font, op["size"] * k))
            out.append("1 0 0 1 {:.2f} {:.2f} Tm".format(x * k, y * k))
            out.append("{} Tj".format(payload))
            out.append("ET")
        elif kind == "rect":
            x = op["x"] * k
            y = (height - op["y"] - op["h"]) * k
            out.append("{:.2f} {:.2f} {:.2f} {:.2f} re".format(x, y, op["w"] * k, op["h"] * k))
            has_fill = op["fill"] is not None
            has_stroke = op["stroke"] is not None
            if has_fill:
                out.append("{} rg".format(_col(op["fill"])))
            if has_stroke:
                out.append("{} RG".format(_col(op["stroke"])))
                out.append("{:.2f} w".format(op["stroke_width"] * k))
            if has_fill and has_stroke:
                out.append("B")
            elif has_fill:
                out.append("f")
            elif has_stroke:
                out.append("S")
        elif kind == "polyline":
            out.append("{} RG".format(_col(op["stroke"])))
            out.append("{:.2f} w".format(op["stroke_width"] * k))
            if op.get("dash"):
                out.append("[{:.2f} {:.2f}] 0 d".format(op["dash"][0] * k, op["dash"][1] * k))
            first = True
            for px, py in op["points"]:
                operator = "m" if first else "l"
                first = False
                out.append("{:.2f} {:.2f} {}".format(px * k, (height - py) * k, operator))
            out.append("S")
            if op.get("dash"):
                out.append("[] 0 d")
        elif kind == "polygon":
            first = True
            for px, py in op["points"]:
                operator = "m" if first else "l"
                first = False
                out.append("{:.2f} {:.2f} {}".format(px * k, (height - py) * k, operator))
            out.append("h")
            has_fill = op["fill"] is not None
            has_stroke = op["stroke"] is not None
            if has_fill:
                out.append("{} rg".format(_col(op["fill"])))
            if has_stroke:
                out.append("{} RG".format(_col(op["stroke"])))
                out.append("{:.2f} w".format(op["stroke_width"] * k))
            if has_fill and has_stroke:
                out.append("B")
            elif has_fill:
                out.append("f")
            elif has_stroke:
                out.append("S")
        elif kind == "circle":
            cx = op["cx"] * k
            cy = (height - op["cy"]) * k
            r = op["r"] * k
            kr = 0.5523 * r
            has_fill = op["fill"] is not None
            has_stroke = op["stroke"] is not None
            if has_fill:
                out.append("{} rg".format(_col(op["fill"])))
            if has_stroke:
                out.append("{} RG".format(_col(op["stroke"])))
                out.append("{:.2f} w".format(op["stroke_width"] * k))
            out.append("{:.2f} {:.2f} m".format(cx - r, cy))
            out.append(
                "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} c".format(
                    cx - r, cy + kr, cx - kr, cy + r, cx, cy + r
                )
            )
            out.append(
                "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} c".format(
                    cx + kr, cy + r, cx + r, cy + kr, cx + r, cy
                )
            )
            out.append(
                "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} c".format(
                    cx + r, cy - kr, cx + kr, cy - r, cx, cy - r
                )
            )
            out.append(
                "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} c".format(
                    cx - kr, cy - r, cx - r, cy - kr, cx - r, cy
                )
            )
            if has_fill and has_stroke:
                out.append("B")
            elif has_fill:
                out.append("f")
            elif has_stroke:
                out.append("S")
        elif kind == "gradient_rect":
            steps = max(2, op["steps"])
            slice_w = op["w"] / steps
            y = (height - op["y"] - op["h"]) * k
            for idx in range(steps):
                frac = float(idx) / max(1, steps - 1)
                color = op["color_fn"](frac)
                out.append("{} rg".format(_col(color)))
                out.append(
                    "{:.2f} {:.2f} {:.2f} {:.2f} re f".format(
                        (op["x"] + idx * slice_w) * k, y, slice_w * k + 0.01, op["h"] * k
                    )
                )
            if op["stroke"]:
                out.append("{} RG".format(_col(op["stroke"])))
                out.append("0.50 w")
                out.append(
                    "{:.2f} {:.2f} {:.2f} {:.2f} re S".format(
                        op["x"] * k, y, op["w"] * k, op["h"] * k
                    )
                )
    out.append("Q")
    return "\n".join(out)


def _tounicode_cmap(entries):
    """ToUnicode CMap stream body: [(gid, char)] -> bfchar blocks of 100."""
    lines = [
        "/CIDInit /ProcSet findresource begin",
        "12 dict begin",
        "begincmap",
        "/CMapName /Adobe-Identity-UCS def",
        "/CMapType 2 def",
        "1 begincodespacerange",
        "<0000> <ffff>",
        "endcodespacerange",
    ]
    for start in range(0, len(entries), 100):
        chunk = entries[start : start + 100]
        lines.append("{} beginbfchar".format(len(chunk)))
        for gid, ch in chunk:
            lines.append("<{:04X}> <{}>".format(gid, ch.encode("utf-16-be").hex().upper()))
        lines.append("endbfchar")
    lines.extend(["endcmap", "CMapName currentdict /CMap defineresource pop", "end", "end"])
    return "\n".join(lines)


def _embedded_cjk_font_objects(used_chars):
    """Build PDF objects for an embedded subset CJK font.

    Returns (dict of objects, next_free_number). Object layout:
      7 Type0 (Identity-H) - 8 CIDFontType2 - 9 FontDescriptor
      10 FontFile2 stream  - 11 ToUnicode stream
    """
    from fontTools.ttLib import TTFont

    ttf_bytes, metrics = _subset_cjk_font(used_chars)
    font = TTFont(metrics["path"], fontNumber=0)
    head = font["head"]
    hhea = font["hhea"]
    upem = float(head.unitsPerEm)
    scale = 1000.0 / upem
    bbox = [
        int(round(head.xMin * scale)),
        int(round(head.yMin * scale)),
        int(round(head.xMax * scale)),
        int(round(head.yMax * scale)),
    ]
    ascent = int(round(hhea.ascent * scale))
    descent = int(round(hhea.descent * scale))

    # Widths and ToUnicode entries keyed by glyph id (retain_gids keeps the
    # numbering of the source font, matching the Identity-H text codes).
    width_entries = []
    tounicode_entries = []
    for ch in used_chars:
        code = ord(ch)
        if code not in metrics["gid"]:
            continue
        gid = metrics["gid"][code]
        advance = int(round(metrics["adv"][code] * scale))
        width_entries.append((gid, advance))
        tounicode_entries.append((gid, ch))
    width_entries.sort()

    w_array = " ".join("{} [{}]".format(gid, adv) for gid, adv in width_entries)
    objects = {}
    objects[7] = (
        "<< /Type /Font /Subtype /Type0 /BaseFont /PVSSubset "
        "/Encoding /Identity-H /DescendantFonts [8 0 R] /ToUnicode 11 0 R >>"
    ).encode("latin-1")
    objects[8] = (
        "<< /Type /Font /Subtype /CIDFontType2 /BaseFont /PVSSubset "
        "/CIDSystemInfo << /Registry (Adobe) /Ordering (Identity) /Supplement 0 >> "
        "/FontDescriptor 9 0 R /DW 1000 /W [{}] /CIDToGIDMap /Identity >>"
    ).format(w_array).encode("latin-1")
    objects[9] = (
        "<< /Type /FontDescriptor /FontName /PVSSubset /Flags 4 "
        "/FontBBox [{0} {1} {2} {3}] /ItalicAngle 0 /Ascent {4} /Descent {5} "
        "/CapHeight {4} /StemV 80 /FontFile2 10 0 R >>"
    ).format(bbox[0], bbox[1], bbox[2], bbox[3], ascent, descent).encode("latin-1")
    objects[10] = (
        b"<< /Length " + str(len(ttf_bytes)).encode("latin-1")
        + b" /Length1 " + str(len(ttf_bytes)).encode("latin-1")
        + b" >>\nstream\n" + ttf_bytes + b"\nendstream"
    )
    tounicode = _tounicode_cmap(tounicode_entries).encode("latin-1")
    objects[11] = (
        b"<< /Length " + str(len(tounicode)).encode("latin-1")
        + b" >>\nstream\n" + tounicode + b"\nendstream"
    )
    return objects, 12


def _fallback_cjk_font_objects():
    """Non-embedded Type0 CJK font (STSong-Light) for fontless environments."""
    objects = {}
    objects[7] = (
        b"<< /Type /Font /Subtype /Type0 /BaseFont /STSong-Light "
        b"/Encoding /UniGB-UCS2-H /DescendantFonts [8 0 R] >>"
    )
    objects[8] = (
        b"<< /Type /Font /Subtype /CIDFontType0 /BaseFont /STSong-Light "
        b"/CIDSystemInfo << /Registry (Adobe) /Ordering (GB1) /Supplement 5 >> /DW 1000 >>"
    )
    return objects, 9


def write_pdf_file(fig, path):
    """Write the figure as a one-page vector PDF.

    Helvetica core fonts carry latin text; CJK text uses an embedded
    subset of an open-licence font when available (any viewer shows it
    correctly), otherwise a non-embedded Type0 reference.
    """
    k = PT_PER_PX
    page_w = max(1.0, fig.width * k)
    page_h = max(1.0, fig.height * k)

    cjk_strings = [
        op["s"] for op in fig.ops if op["op"] == "text" and _has_nonlatin(op["s"])
    ]
    if not cjk_strings:
        cjk_mode = "none"
    elif _cjk_metrics():
        cjk_mode = "gid"
    else:
        cjk_mode = "ucs2"

    content = render_pdf_stream(fig, cjk_mode=cjk_mode).encode("latin-1", "replace")
    objects = {}
    objects[1] = b"<< /Type /Catalog /Pages 2 0 R >>"
    objects[2] = b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>"
    resources = "<< /Font << /F1 5 0 R /F2 6 0 R"
    if cjk_mode != "none":
        resources += " /F3 7 0 R"
    resources += " >> >>"
    objects[3] = (
        "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {w:.2f} {h:.2f}] "
        "/Resources {r} /Contents 4 0 R >>"
    ).format(w=page_w, h=page_h, r=resources).encode("latin-1")
    objects[4] = (
        b"<< /Length " + str(len(content)).encode("latin-1") + b" >>\nstream\n" + content + b"\nendstream"
    )
    objects[5] = b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>"
    objects[6] = b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>"
    if cjk_mode == "gid":
        cjk_objects, count = _embedded_cjk_font_objects(set("".join(cjk_strings)))
    elif cjk_mode == "ucs2":
        cjk_objects, count = _fallback_cjk_font_objects()
    else:
        cjk_objects, count = {}, 7
    objects.update(cjk_objects)

    buf = bytearray()
    buf += b"%PDF-1.4\n"
    offsets = {}
    for num in sorted(objects):
        offsets[num] = len(buf)
        buf += "{} 0 obj\n".format(num).encode("latin-1") + objects[num] + b"\nendobj\n"
    xref_pos = len(buf)
    buf += "xref\n0 {}\n".format(count).encode("latin-1")
    buf += b"0000000000 65535 f \n"
    for num in sorted(objects):
        buf += "{:010d} 00000 n \n".format(offsets[num]).encode("latin-1")
    buf += "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{}\n%%EOF\n".format(
        count, xref_pos
    ).encode("latin-1")
    path = Path(path)
    path.write_bytes(bytes(buf))
    return path


# ---------------------------------------------------------------------------
# Output bundle (interface kept from the per-script implementations)
# ---------------------------------------------------------------------------


def write_output_bundle(out_path, html_doc, write_svg=True, write_pdf=True):
    """Write the figure bundle next to ``out_path``.

    Writes ``out_path`` (.html) and, unless ``write_pdf`` is False, a native
    vector ``.pdf`` rendered from the same Figure. ``write_svg`` is accepted
    for backward compatibility but ignored: standalone SVG output was
    replaced by the PDF export. Returns ``(html_path, None, pdf_path)``.
    """
    out_path = Path(out_path)
    out_path.write_text(html_doc, encoding="utf-8")
    pdf_path = None
    if write_pdf:
        figure = getattr(html_doc, "figure", None)
        if figure is None:
            raise RuntimeError(
                "PDF export needs an HtmlDoc produced by render_html (native vector export)."
            )
        pdf_path = out_path.with_suffix(".pdf")
        write_pdf_file(figure, pdf_path)
    return out_path, None, pdf_path
