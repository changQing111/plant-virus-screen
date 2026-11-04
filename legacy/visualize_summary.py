#!/usr/bin/env python3
"""Report-level summary dashboard for plant-virus-screen.

One three-panel figure rendered through the shared layer (self-contained
HTML plus native vector PDF):

1. Evidence grade distribution (A/B/C/D, log scale, counts and shares)
2. Top reference-related groups by contig count, stacked by grade
   composition
3. Coverage overview scatter: MAPQ>=20 mean depth vs breadth>=1x for
   A/B/C candidates, with the report mapping thresholds marked

Reads the TSVs that report.py already wrote into the report directory,
so it can be regenerated standalone at any time:

    python3 visualize_summary.py --report-dir <report_dir>
"""

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

_LEGACY_DIR = str(Path(__file__).resolve().parent)
if _LEGACY_DIR not in sys.path:
    sys.path.insert(0, _LEGACY_DIR)

import render_common as rc

GRADE_ORDER = ["A_supported", "B_supported", "C_provisional", "D_conflicting"]
GRADE_SHORT = {
    "A_supported": "A",
    "B_supported": "B",
    "C_provisional": "C",
    "D_conflicting": "D",
}
GRADE_COLOR = {
    "A_supported": "#0F6E56",
    "B_supported": "#5DCAA5",
    "C_provisional": "#EF9F27",
    "D_conflicting": "#B4B2A9",
}


def _log10(value):
    if value <= 0:
        return 0.0
    return math.log10(value)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build the report summary dashboard (HTML + native vector PDF) from report TSVs."
    )
    parser.add_argument("--report-dir", required=True, help="Report directory with annotated_candidates.tsv")
    parser.add_argument("--breadth", type=float, default=0.8, help="MAPQ>=20 breadth>=1x threshold line")
    parser.add_argument("--depth", type=float, default=5.0, help="MAPQ>=20 mean depth threshold line")
    parser.add_argument("--no-pdf", action="store_true", help="Only write HTML.")
    return parser.parse_args()


def _read_tsv(path):
    rows = []
    if not path or not Path(path).is_file():
        return rows
    with Path(path).open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(row)
    return rows


def _short_contig_label(contig_id):
    match = re.search(r"NODE_(\d+)_length_(\d+)", contig_id)
    if match:
        length = int(match.group(2))
        if length >= 1000:
            return "NODE_{}/{:.1f}kb".format(match.group(1), length / 1000.0)
        return "NODE_{}/{}bp".format(match.group(1), length)
    return contig_id[:16]


def _truncate(text, limit):
    if len(text) <= limit:
        return text
    return text[: limit - 3] + "..."


def _fmt_count(value):
    if value >= 10000:
        return "{:.1f}k".format(value / 1000.0)
    if value >= 1000:
        return "{:.2f}k".format(value / 1000.0)
    return str(value)


def _panel_grade_distribution(fig, counts, top_y, base_y, x_left, x_right):
    label_style = {"size": 12, "fill": "#222222"}
    axis_style = {"size": 11, "fill": "#555555"}
    small_style = {"size": 10, "fill": "#555555"}
    total = sum(counts.get(g, 0) for g in GRADE_ORDER) or 1

    fig.text(24, top_y - 10, "Evidence grade distribution (log scale)", bold=True, **{"size": 13, "fill": "#222222"})

    max_count = max([counts.get(g, 0) for g in GRADE_ORDER] + [1])
    log_max = 1.0
    while 10 ** log_max <= max_count:
        log_max += 1
    plot_h = base_y - (top_y + 16)
    slot = (x_right - x_left) / float(len(GRADE_ORDER))
    bar_w = min(140.0, slot * 0.5)

    # Horizontal gridlines at powers of ten.
    exponent = 0
    while 10 ** exponent <= 10 ** (log_max - 1):
        y = base_y - (exponent / log_max) * plot_h
        fig.polyline([(x_left, y), (x_right, y)], stroke="#e0e0e0", stroke_width=0.5)
        label = 1 if exponent == 0 else _fmt_count(10 ** exponent)
        fig.text(x_left - 14, y + 4, label, anchor="end", **small_style)
        exponent += 1
    fig.polyline([(x_left, base_y), (x_right, base_y)], stroke="#555555", stroke_width=1.0)

    for index, grade in enumerate(GRADE_ORDER):
        count = counts.get(grade, 0)
        center = x_left + slot * (index + 0.5)
        if count <= 0:
            bar_h = 2.0
        elif count == 1:
            bar_h = 6.0
        else:
            bar_h = max(2.0, (_log10(count) / log_max) * plot_h)
        fig.rect(center - bar_w / 2, base_y - bar_h, bar_w, bar_h, fill=GRADE_COLOR[grade])
        fig.text(center, base_y - bar_h - 8, str(count), anchor="middle", **label_style)
        share = 100.0 * count / total
        fig.text(center, base_y - bar_h - 24, "({:.2f}%)".format(share), anchor="middle", **small_style)
        fig.text(center, base_y + 18, GRADE_SHORT[grade], anchor="middle", **label_style)
        fig.text(center, base_y + 34, _grade_hint(grade), anchor="middle", **small_style)


def _grade_hint(grade):
    return {
        "A_supported": "multi-evidence + mapping",
        "B_supported": "multi-evidence",
        "C_provisional": "single/unresolved",
        "D_conflicting": "cellular conflict",
    }.get(grade, "")


def _panel_group_stack(fig, groups, top_y, max_rows, x_left, x_right):
    label_style = {"size": 11, "fill": "#222222"}
    small_style = {"size": 10, "fill": "#555555"}

    fig.text(24, top_y - 10, "Top reference-related groups by contig count (grade composition)", bold=True, **{"size": 13, "fill": "#222222"})

    ordered = sorted(groups, key=lambda g: int(float(g.get("contig_count") or 0)), reverse=True)
    ordered = ordered[:max_rows]
    if not ordered:
        fig.text(24, top_y + 20, "No reference-related groups available.", **small_style)
        return top_y + 40

    bar_x0 = x_left + 40
    bar_x1 = x_right
    max_count = max(int(float(g["contig_count"])) for g in ordered) or 1
    row_h = 24
    legend_y = top_y + 12
    for index, grade in enumerate(GRADE_ORDER):
        lx = bar_x1 - 460 + index * 115
        fig.rect(lx, legend_y - 8, 10, 10, fill=GRADE_COLOR[grade])
        fig.text(lx + 14, legend_y, GRADE_SHORT[grade], **small_style)

    y = legend_y + 18
    for group in ordered:
        name = _truncate((group.get("reference_related_group") or "").strip(), 52)
        count = int(float(group["contig_count"]))
        grade_counts = {}
        try:
            grade_counts = json.loads(group.get("grade_counts") or "{}")
        except ValueError:
            grade_counts = {}
        fig.text(x_left, y + 6, name, **label_style)
        total_w = (bar_x1 - bar_x0) * (count / float(max_count))
        x_cursor = bar_x0
        for grade in GRADE_ORDER:
            part = int(grade_counts.get(grade, 0) or 0)
            if part <= 0:
                continue
            part_w = total_w * (part / float(count))
            fig.rect(x_cursor, y - 6, part_w, 12, fill=GRADE_COLOR[grade])
            x_cursor += part_w
        fig.text(min(bar_x1 + 8, x_cursor + 6), y + 4, str(count), **small_style)
        y += row_h
    return y


def _panel_coverage_scatter(fig, rows, top_y, base_y, x_left, x_right, breadth_threshold, depth_threshold):
    label_style = {"size": 12, "fill": "#222222"}
    axis_style = {"size": 11, "fill": "#555555"}
    small_style = {"size": 10, "fill": "#555555"}

    fig.text(24, top_y - 10, "Coverage overview (MAPQ>=20): mean depth vs breadth>=1x, A/B/C candidates", bold=True, **{"size": 13, "fill": "#222222"})

    points = []
    for row in rows:
        grade = (row.get("evidence_grade") or "").strip()
        if grade not in ("A_supported", "B_supported", "C_provisional"):
            continue
        try:
            mean_depth = float(row.get("mapq20_mean_depth") or "")
            breadth = float(row.get("mapq20_breadth_1x") or "")
        except ValueError:
            continue
        points.append((row.get("contig_id") or "", grade, mean_depth, breadth))

    if not points:
        fig.text(24, top_y + 20, "No MAPQ>=20 coverage data for A/B/C candidates; panel skipped.", **small_style)
        return

    plot_h = base_y - (top_y + 16)
    x0 = x_left + 40
    x1 = x_right
    max_depth = max([p[2] for p in points] + [depth_threshold, 1.0])
    log_x_max = _log10(max_depth + 1.0) or 1.0

    def x_of(depth):
        return x0 + (_log10(depth + 1.0) / log_x_max) * (x1 - x0)

    def y_of(breadth):
        return base_y - max(0.0, min(1.0, breadth)) * plot_h

    fig.rect(x0, top_y + 16, x1 - x0, plot_h, fill="#fcfcfc", stroke="#d0d0d0", stroke_width=1.0)

    # Threshold lines (report mapping criteria).
    if 0.0 < breadth_threshold < 1.0:
        y_th = y_of(breadth_threshold)
        fig.polyline([(x0, y_th), (x1, y_th)], stroke="#d85a30", stroke_width=1.0, dash=(6, 4))
        fig.text(x1 + 6, y_th + 4, "breadth {:.2f}".format(breadth_threshold), **small_style)
    x_th = x_of(depth_threshold)
    if x0 < x_th < x1:
        fig.polyline([(x_th, top_y + 16), (x_th, base_y)], stroke="#d85a30", stroke_width=1.0, dash=(6, 4))
        fig.text(x_th + 4, top_y + 30, "depth {}".format(_fmt_count(depth_threshold)), **small_style)

    # Axis ticks: x powers of ten on the log scale, y quarters.
    exponent = 0
    while 10 ** exponent <= max_depth:
        tick_x = x_of(10 ** exponent - 1.0)
        if tick_x <= x1:
            fig.text(tick_x, base_y + 18, _fmt_count(10 ** exponent), anchor="middle", **axis_style)
        exponent += 1
    fig.text((x0 + x1) / 2.0, base_y + 38, "mean depth (MAPQ>=20, log scale)", anchor="middle", **axis_style)
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        fig.text(x0 - 10, y_of(frac) + 4, "{:.2f}".format(frac), anchor="end", **axis_style)
    fig.text(x0 - 56, top_y + 26, "breadth >=1x", **axis_style)

    for contig, grade, mean_depth, breadth in points:
        color = GRADE_COLOR[grade]
        px = x_of(mean_depth)
        py = y_of(breadth)
        fig.circle(px, py, 5, fill=color, title="{} | {} | depth {:.1f} | breadth {:.3f}".format(contig, grade, mean_depth, breadth))
        if grade == "A_supported":
            fig.text(px + 8, py + 4, _short_contig_label(contig), **{"size": 9, "fill": "#333333"})

    # Legend inside the plot, top-right.
    legend_x = x1 - 300
    legend_y = top_y + 40
    for index, grade in enumerate(["A_supported", "B_supported", "C_provisional"]):
        fig.circle(legend_x + index * 95, legend_y, 5, fill=GRADE_COLOR[grade])
        fig.text(legend_x + index * 95 + 10, legend_y + 4, GRADE_SHORT[grade], **small_style)
    fig.text(x0 + 8, top_y + 16 + plot_h - 10, "Dashed lines mark the report mapping thresholds; grade A additionally requires multi-source viral evidence.", **small_style)


def _panel_glossary(fig, top_y, breadth_threshold, depth_threshold):
    """Bilingual glossary block at the bottom of the dashboard."""
    fig.text(24, top_y, "Glossary 术语表", bold=True, **{"size": 13, "fill": "#222222"})
    lines = [
        "Grade A | A 级：至少两类证据（nt / 蛋白 / RdRP 结构域）、竞争最佳命中为病毒，且回贴广度≥1x 不低于 {:.2f}、平均深度不低于 {}".format(breadth_threshold, _fmt_count(depth_threshold)),
        "Grade B | B 级：同样多类证据且竞争最佳命中为病毒，但回贴验证缺失或未达标",
        "Grade C | C 级：仅有单一类证据，或竞争情况未解决",
        "Grade D | D 级：竞争最佳命中为细胞来源（宿主 / 细菌 / 真菌），大概率非病毒序列",
        "breadth >=1x | 广度：contig 上被至少 1 条高质量回贴 reads 覆盖的位置比例（面板三纵轴）",
        "mean depth | 平均深度：逐位置回贴深度的平均值（面板三横轴，对数刻度）",
        "MAPQ>=20 | 仅统计比对质量不低于 20 的回贴 reads，排除错配与多位点比对的噪声",
        "Reference-related group | 参考相关分组：推断与同一条参考病毒相关的全部 contig（面板二每行一组）",
        "Log scale | 对数轴：面板一纵轴与面板三横轴按 10 的幂刻度，用于同时展示悬殊的数量级",
        "Grades are heuristic review priorities, not infection or species confirmation | 分级仅是启发式复核优先级，不构成感染、物种或新病毒的确认",
    ]
    y = top_y + 26
    for line in lines:
        fig.text(24, y, line, **{"size": 11, "fill": "#555555"})
        y += 18
    return y + 16


def build_figure(rows, groups, breadth_threshold, depth_threshold):
    counts = {}
    for row in rows:
        grade = (row.get("evidence_grade") or "").strip() or "unknown"
        counts[grade] = counts.get(grade, 0) + 1
    total = len(rows)

    fig = rc.Figure(1400, 1120)
    fig.text(24, 34, "Report summary", bold=True, **{"size": 15, "fill": "#222222"})
    meta_bits = ["Candidates: {}".format(total)]
    meta_bits.extend("{}: {}".format(GRADE_SHORT.get(g, g), counts.get(g, 0)) for g in GRADE_ORDER)
    meta_bits.append("Groups: {}".format(len(groups)))
    meta_text = " | ".join(meta_bits)
    fig.text(24, 56, meta_text, **{"size": 11, "fill": "#666666"})

    _panel_grade_distribution(fig, counts, 90, 270, 120, 1330)
    next_y = _panel_group_stack(fig, groups, 320, 15, 24, 1330)
    scatter_top = max(760, next_y + 40)
    _panel_coverage_scatter(fig, rows, scatter_top, scatter_top + 260, 24, 1330, breadth_threshold, depth_threshold)
    glossary_bottom = _panel_glossary(fig, scatter_top + 260 + 70, breadth_threshold, depth_threshold)
    fig.height = max(1120, glossary_bottom)
    return fig, meta_text


def build_summary(report_dir, breadth_threshold=0.8, depth_threshold=5.0, write_pdf=True):
    """Build report_summary.html (+PDF) from the report TSVs.

    Returns (html_path, pdf_path) or None when no annotated candidates
    table exists.
    """
    report_dir = Path(report_dir)
    rows = _read_tsv(report_dir / "annotated_candidates.tsv")
    if not rows:
        return None
    groups = _read_tsv(report_dir / "virus_evidence_summary.tsv")
    fig, meta_text = build_figure(rows, groups, breadth_threshold, depth_threshold)
    doc = rc.render_html(fig, "Report summary", "Report summary", meta_text)
    return rc.write_output_bundle(report_dir / "report_summary.html", doc, write_pdf=write_pdf)


def main():
    args = parse_args()
    result = build_summary(
        args.report_dir, breadth_threshold=args.breadth, depth_threshold=args.depth, write_pdf=not args.no_pdf
    )
    if result is None:
        raise SystemExit("No annotated_candidates.tsv under {}; nothing to summarize.".format(args.report_dir))
    html_path, svg_path, pdf_path = result
    print("Wrote HTML:", html_path.resolve())
    if pdf_path is not None:
        print("Wrote PDF:", pdf_path.resolve())


if __name__ == "__main__":
    main()
