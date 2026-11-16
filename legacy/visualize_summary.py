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


def _wrap_text(text, width, size=11):
    """Wrap by rendered width, including unbroken IDs; retain every character."""
    lines = []
    remaining = str(text).strip()
    while remaining:
        if rc.text_width(remaining, size) <= width:
            lines.append(remaining)
            break
        cut = 1
        while cut < len(remaining) and rc.text_width(remaining[:cut + 1], size) <= width:
            cut += 1
        space = remaining.rfind(" ", 0, cut + 1)
        if space > cut // 2:
            cut = space
        lines.append(remaining[:cut].rstrip())
        remaining = remaining[cut:].lstrip()
    return lines or [""]


def _panel_group_stack(fig, groups, top_y, max_rows, x_left, x_right):
    small_style = {"size": 10, "fill": "#555555"}
    fig.text(24, top_y - 10, "Top reference-related groups by contig count (grade composition)",
             size=13, fill="#222222", bold=True)
    ordered = sorted(groups, key=lambda g: int(float(g.get("contig_count") or 0)), reverse=True)[:max_rows]
    if not ordered:
        fig.text(24, top_y + 20, "No reference-related groups available.", **small_style)
        return top_y + 40

    # Independent label, bar and count columns; widths use the rendering metrics.
    bar_x0 = max(x_left + 410, 440)
    bar_x1 = x_right - 55
    max_count = max(int(float(g["contig_count"])) for g in ordered) or 1
    legend_y = top_y + 14
    for index, grade in enumerate(GRADE_ORDER):
        lx = bar_x0 + index * 90
        fig.rect(lx, legend_y - 8, 10, 10, fill=GRADE_COLOR[grade])
        fig.text(lx + 15, legend_y, GRADE_SHORT[grade], **small_style)
    fig.text(bar_x1 + 10, legend_y, "Contigs", **small_style)
    y = legend_y + 26
    for group in ordered:
        name = (group.get("reference_related_group") or "").strip()
        lines = _wrap_text(name, bar_x0 - x_left - 24)
        row_h = max(30, len(lines) * 15 + 12)
        center = y + row_h / 2
        first_baseline = center - (len(lines) - 1) * 7.5 + 4
        for i, line in enumerate(lines):
            fig.text(x_left, first_baseline + 15 * i, line, size=11, fill="#222222", title=name)
        count = int(float(group["contig_count"]))
        try:
            grade_counts = json.loads(group.get("grade_counts") or "{}")
        except ValueError:
            grade_counts = {}
        total_w = (bar_x1 - bar_x0) * (count / float(max_count))
        x_cursor = bar_x0
        for grade in GRADE_ORDER:
            part = int(grade_counts.get(grade, 0) or 0)
            if part <= 0 or count <= 0:
                continue
            part_w = total_w * part / count
            fig.rect(x_cursor, center - 6, part_w, 12, fill=GRADE_COLOR[grade],
                     title="{} | {}: {}".format(name, GRADE_SHORT[grade], part))
            x_cursor += part_w
        fig.text(bar_x1 + 10, center + 4, str(count), **small_style)
        y += row_h
    return y


def _panel_coverage_scatter(fig, rows, top_y, base_y, x_left, x_right, breadth_threshold, depth_threshold):
    axis_style = {"size": 11, "fill": "#555555"}
    small_style = {"size": 10, "fill": "#555555"}
    fig.text(24, top_y - 10, "Coverage overview (MAPQ>=20): mean depth vs breadth>=1x, A/B/C candidates",
             size=13, fill="#222222", bold=True)
    points = []
    for row in rows:
        grade = (row.get("evidence_grade") or "").strip()
        if grade not in GRADE_ORDER[:3]:
            continue
        try:
            mean_depth = float(row.get("mapq20_mean_depth") or "")
            breadth = float(row.get("mapq20_breadth_1x") or "")
        except ValueError:
            continue
        if not (math.isfinite(mean_depth) and math.isfinite(breadth) and mean_depth >= 0 and 0 <= breadth <= 1):
            continue
        points.append((row.get("contig_id") or "", grade, mean_depth, breadth))
    if not points:
        fig.text(24, top_y + 20, "No valid MAPQ>=20 coverage data for A/B/C candidates; panel skipped.", **small_style)
        return top_y + 50

    x0, x1 = max(x_left + 80, 105), x_right - 40
    plot_top = top_y + 70
    plot_h = base_y - plot_top
    max_depth = max([p[2] for p in points] + [depth_threshold, 1.0])
    # Add headroom so breadth=1 and maximum-depth points remain inside the axes.
    log_x_max = _log10(max_depth + 1.0) * 1.04 or 1.0
    def x_of(depth):
        return x0 + (_log10(depth + 1.0) / log_x_max) * (x1 - x0)
    def y_of(breadth):
        return base_y - (breadth + 0.025) / 1.05 * plot_h

    legend_y = top_y + 15
    for index, grade in enumerate(GRADE_ORDER[:3]):
        lx = x0 + index * 90
        fig.circle(lx, legend_y - 3, 5, fill=GRADE_COLOR[grade])
        fig.text(lx + 12, legend_y, GRADE_SHORT[grade], **small_style)
    fig.text(x0 + 310, legend_y,
             "Dashed thresholds: breadth {:.2f}; mean depth {}".format(breadth_threshold, _fmt_count(depth_threshold)), **small_style)
    fig.text(x0, top_y + 47, "Breadth >=1x (MAPQ>=20)", **axis_style)
    fig.rect(x0, plot_top, x1 - x0, plot_h, fill="#fcfcfc", stroke="#d0d0d0", stroke_width=1)
    if 0 <= breadth_threshold <= 1:
        y_th = y_of(breadth_threshold)
        fig.polyline([(x0, y_th), (x1, y_th)], stroke="#d85a30", stroke_width=1, dash=(6, 4))
    x_th = x_of(depth_threshold)
    if x0 < x_th < x1:
        fig.polyline([(x_th, plot_top), (x_th, base_y)], stroke="#d85a30", stroke_width=1, dash=(6, 4))
    # Coordinates use log10(depth+1); labels show actual depth, including zero.
    ticks = [0.0]
    exponent = 0
    while 10 ** exponent <= max_depth:
        ticks.append(float(10 ** exponent)); exponent += 1
    last_right = -math.inf
    for depth in ticks:
        px = x_of(depth); label = _fmt_count(int(depth))
        half = rc.text_width(label, 11) / 2
        if px - half > last_right + 12:
            fig.text(px, base_y + 22, label, anchor="middle", **axis_style)
            last_right = px + half
    fig.text((x0 + x1) / 2, base_y + 44, "Mean depth (MAPQ>=20; log10(depth + 1) spacing)", anchor="middle", **axis_style)
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        fig.text(x0 - 12, y_of(frac) + 4, "{:.2f}".format(frac), anchor="end", **axis_style)

    # Dense scatter: draw C then B then A; full IDs remain in SVG tooltips.
    for contig, grade, mean_depth, breadth in sorted(points, key=lambda p: GRADE_ORDER.index(p[1]), reverse=True):
        fig.circle(x_of(mean_depth), y_of(breadth), 4, fill=GRADE_COLOR[grade],
                   title="{} | {} | depth {:.2f} | breadth {:.4f}".format(contig, grade, mean_depth, breadth))
    y = base_y + 70
    fig.text(24, y, "All A/B/C points are retained. Hover over points in HTML for IDs and metrics; grade A candidates are listed below.", **small_style)
    y += 26
    supported = [p for p in points if p[1] == "A_supported"]
    if supported:
        fig.text(24, y, "Grade A candidate details (not point annotations)", size=11, bold=True)
        y += 22
        fig.text(24, y, "Contig ID", **small_style)
        fig.text(1100, y, "Mean depth", anchor="end", **small_style)
        fig.text(1300, y, "Breadth >=1x", anchor="end", **small_style)
        y += 20
        for contig, grade, mean_depth, breadth in sorted(supported):
            lines = _wrap_text(contig, 970, 10)
            for i, line in enumerate(lines):
                fig.text(24, y + i * 14, line, **small_style)
            fig.text(1100, y, "{:.2f}".format(mean_depth), anchor="end", **small_style)
            fig.text(1300, y, "{:.4f}".format(breadth), anchor="end", **small_style)
            y += max(20, len(lines) * 14 + 6)
    return y + 16


def _panel_glossary(fig, top_y, breadth_threshold, depth_threshold):
    fig.text(24, top_y, "Glossary 术语表", bold=True, size=13, fill="#222222")
    lines = [
        "Grade A | A 级：至少两类证据且包含 nt、竞争汇总为 viral，MAPQ20 广度不低于 {:.2f}、平均深度不低于 {}".format(breadth_threshold, _fmt_count(depth_threshold)),
        "Grade B | B 级：至少两类证据，且不满足 A 或 D；竞争未分类或缺少 nt 的候选也可为 B",
        "Grade C | C 级：其余候选，通常为单类证据或未解决证据",
        "Grade D | D 级：竞争汇总为 cellular，表示来源冲突，不能据此排除病毒",
        "Discovery priority | 新病毒复核请同时查看 discovery_priority；mixed / ambiguous 对应 review_conflict，不自动归为 D",
        "breadth >=1x | 广度：contig 中至少被 1 条 MAPQ>=20 回贴 read 覆盖的位置比例",
        "mean depth | 平均深度：完整 contig 上逐位置深度的平均值；横轴位置为 log10(depth + 1)，刻度显示实际深度",
        "MAPQ>=20 | 回贴比对质量过滤，不保证排除所有错配、多位点比对或污染",
        "Reference-related group | 参考相关分组：按推断参考名称汇总的 contig；不是已确认病毒物种数",
        "Grades and priorities are heuristic; they do not confirm infection, species or novelty | 分级及优先级不能确认感染、物种或新病毒",
    ]
    y = top_y + 26
    for line in lines:
        for wrapped in _wrap_text(line, fig.width - 48, 11):
            fig.text(24, y, wrapped, size=11, fill="#555555")
            y += 18
    return y + 16


def build_figure(rows, groups, breadth_threshold, depth_threshold):
    counts = {}
    for row in rows:
        grade = (row.get("evidence_grade") or "").strip() or "unknown"
        counts[grade] = counts.get(grade, 0) + 1
    fig = rc.Figure(1400, 1120)
    fig.text(24, 34, "Report summary", bold=True, size=15, fill="#222222")
    meta_bits = ["Candidates: {}".format(len(rows))]
    meta_bits.extend("{}: {}".format(GRADE_SHORT.get(g, g), counts.get(g, 0)) for g in GRADE_ORDER)
    meta_bits.append("Groups: {}".format(len(groups)))
    meta_text = " | ".join(meta_bits)
    fig.text(24, 56, meta_text, size=11, fill="#666666")
    _panel_grade_distribution(fig, counts, 90, 270, 120, 1330)
    next_y = _panel_group_stack(fig, groups, 355, 15, 24, 1330)
    scatter_top = next_y + 60
    scatter_bottom = _panel_coverage_scatter(fig, rows, scatter_top, scatter_top + 360,
                                              24, 1330, breadth_threshold, depth_threshold)
    glossary_bottom = _panel_glossary(fig, scatter_bottom + 24, breadth_threshold, depth_threshold)
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
