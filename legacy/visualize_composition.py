#!/usr/bin/env python3
"""VMR-based composition charts for plant-virus-screen reports.

One figure (self-contained HTML + native vector PDF) with three stacked
composition panels restricted to the review-relevant A/B/C candidates:

1. Host source composition (VMR reference host of the best hit)
2. Genome composition / Baltimore group (VMR genome type of the best hit)
3. Plant-virus status (this pipeline's own host-mapping category)

Bars are stacked by evidence grade so composition and evidence quality
are read together. Chinese labels are rendered directly (the PDF backend
emits CJK through a Type0 font); the figure footer carries the VMR
annotation-scope disclaimer.

Reads the TSVs that report.py already wrote, so it can be regenerated
standalone:

    python3 visualize_composition.py --report-dir <report_dir>
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

_LEGACY_DIR = str(Path(__file__).resolve().parent)
if _LEGACY_DIR not in sys.path:
    sys.path.insert(0, _LEGACY_DIR)

import render_common as rc

REVIEW_GRADES = ["A_supported", "B_supported", "C_provisional"]
GRADE_SHORT = {"A_supported": "A", "B_supported": "B", "C_provisional": "C"}
GRADE_COLOR = {"A_supported": "#0F6E56", "B_supported": "#5DCAA5", "C_provisional": "#EF9F27"}
GRADE_LABEL = {
    "A_supported": "A 多证据+回贴达标",
    "B_supported": "B 多证据",
    "C_provisional": "C 单一/未解决",
}

HOST_SOURCE_ZH = {
    "plants": "植物",
    "bacteria": "细菌",
    "fungi": "真菌",
    "protists": "原生生物",
    "vertebrates": "脊椎动物",
    "invertebrates": "无脊椎动物",
    "archaea": "古菌",
    "algae": "藻类",
}

GENOME_ZH = {
    "dsDNA": "双链DNA",
    "ssDNA": "单链DNA",
    "ssDNA(+)": "正单链DNA",
    "ssDNA(-)": "负单链DNA",
    "dsRNA": "双链RNA",
    "ssRNA(+)": "正单链RNA",
    "ssRNA(-)": "负单链RNA",
    "ssRNA-RT": "单链RNA逆转录",
    "dsDNA-RT": "双链DNA逆转录",
}

STATUS_ZH = {
    "plant virus": "植物病毒",
    "plant_host_not_documented_by_this_mapping": "植物宿主未记载",
    "unknown": "未注释",
}

UNANNOTATED = "未注释"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build the A/B/C composition figure (HTML + native vector PDF) from report TSVs."
    )
    parser.add_argument("--report-dir", required=True, help="Report directory with annotated_candidates.tsv")
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


def _grade_counts_per_value(rows, field, translator):
    """{value_label: {grade: count}} over the A/B/C rows of one field."""
    counts = defaultdict(Counter)
    for row in rows:
        if (row.get("evidence_grade") or "").strip() not in REVIEW_GRADES:
            continue
        raw = (row.get(field) or "").strip()
        label = translator(raw) if raw else UNANNOTATED
        counts[label][row["evidence_grade"]] += 1
    return counts


def _genome_label(row):
    composition = (row.get("genome_composition") or "").strip()
    if not composition:
        return UNANNOTATED
    chinese = (row.get("virus_genome_type_cn") or "").strip() or GENOME_ZH.get(composition, "")
    baltimore = (row.get("baltimore_group") or "").strip()
    bits = [composition, chinese]
    if baltimore:
        bits.append("Baltimore " + baltimore)
    return " · ".join(bits)


def _host_label(raw):
    if not raw:
        return UNANNOTATED
    zh = HOST_SOURCE_ZH.get(raw, raw)
    return "{} {}".format(zh, raw)


def _status_label(raw):
    if not raw:
        return UNANNOTATED
    return STATUS_ZH.get(raw, raw)


def _stacked_panel(fig, top_y, title, counts, x_left, x_right):
    """Horizontal stacked bars sorted by total; returns bottom y."""
    label_style = {"size": 11, "fill": "#222222"}
    small_style = {"size": 10, "fill": "#555555"}

    fig.text(24, top_y, title, bold=True, **{"size": 13, "fill": "#222222"})
    if not counts:
        fig.text(24, top_y + 22, "No annotated candidates in this category.", **small_style)
        return top_y + 46

    ordered = sorted(counts.items(), key=lambda kv: sum(kv[1].values()), reverse=True)
    bar_x0 = x_left + 340
    bar_x1 = x_right - 60
    max_total = max(sum(grades.values()) for _, grades in ordered) or 1
    row_h = 30
    bar_h = 14

    y = top_y + 24
    for label, grades in ordered:
        total = sum(grades.values())
        fig.text(x_left, y + 5, label, **label_style)
        bar_w = (bar_x1 - bar_x0) * (total / float(max_total))
        x_cursor = bar_x0
        for grade in REVIEW_GRADES:
            part = grades.get(grade, 0)
            if part <= 0:
                continue
            part_w = bar_w * (part / float(total))
            fig.rect(x_cursor, y - 3, part_w, bar_h, fill=GRADE_COLOR[grade])
            if part_w >= 46:
                fig.text(
                    x_cursor + part_w / 2.0,
                    y + 8,
                    str(part),
                    anchor="middle",
                    **{"size": 9, "fill": "#ffffff" if grade == "A_supported" else "#222222"}
                )
            x_cursor += part_w
        fig.text(min(bar_x1 + 10, x_cursor + 8), y + 5, str(total), **small_style)
        y += row_h
    return y + 16


def build_figure(rows):
    review = [r for r in rows if (r.get("evidence_grade") or "").strip() in REVIEW_GRADES]
    width = 1400

    host_counts = _grade_counts_per_value(review, "host_source", _host_label)
    genome_counts = defaultdict(Counter)
    for row in review:
        label = _genome_label(row)
        genome_counts[label][(row.get("evidence_grade") or "").strip()] += 1
    status_counts = _grade_counts_per_value(review, "plant_virus_status", _status_label)

    row_counts = (
        max(4, len(host_counts)) + max(4, len(genome_counts)) + max(3, len(status_counts)) + 3
    )
    height = 170 + row_counts * 30 + 90
    fig = rc.Figure(width, height)

    fig.text(24, 34, "A/B/C 候选分类组成 Composition of review candidates", bold=True, **{"size": 15, "fill": "#222222"})
    meta_bits = ["A/B/C candidates: {}".format(len(review))]
    meta_bits.extend(
        "{} {}".format(GRADE_SHORT[g], sum(1 for r in review if (r.get("evidence_grade") or "").strip() == g))
        for g in REVIEW_GRADES
    )
    meta_text = " | ".join(meta_bits)
    fig.text(24, 56, meta_text, **{"size": 11, "fill": "#666666"})

    # Grade legend.
    legend_x = 900
    for index, grade in enumerate(REVIEW_GRADES):
        lx = legend_x + index * 160
        fig.rect(lx, 74, 12, 12, fill=GRADE_COLOR[grade])
        fig.text(lx + 18, 82, GRADE_LABEL[grade], **{"size": 10, "fill": "#555555"})

    y = 110
    y = _stacked_panel(fig, y, "宿主来源组成 Host source (VMR reference host of best hit)", host_counts, 24, 1330)
    y = _stacked_panel(fig, y + 14, "基因组组成 Genome composition / Baltimore group", genome_counts, 24, 1330)
    y = _stacked_panel(fig, y + 14, "植物病毒状态 Plant-virus status (pipeline host mapping)", status_counts, 24, 1330)

    fig.text(
        24,
        y + 20,
        "注：宿主与基因组字段仅描述参考病毒（legacy accession/name/genus 映射），不代表候选的实际宿主；"
        "未注释 = 该候选无对应 VMR 映射。分级是复核优先级，不是感染或物种确认。",
        **{"size": 10, "fill": "#777777"}
    )
    fig.height = y + 46
    return fig, meta_text


def build_composition(report_dir, write_pdf=True):
    """Build composition_summary.html (+PDF) from the report TSVs.

    Returns (html_path, pdf_path) or None when there is nothing to show.
    """
    report_dir = Path(report_dir)
    rows = _read_tsv(report_dir / "annotated_candidates.tsv")
    if not rows:
        return None
    review = [r for r in rows if (r.get("evidence_grade") or "").strip() in REVIEW_GRADES]
    if not review:
        return None
    fig, meta_text = build_figure(rows)
    doc = rc.render_html(fig, "Composition of review candidates", "A/B/C 候选分类组成", meta_text)
    return rc.write_output_bundle(report_dir / "composition_summary.html", doc, write_pdf=write_pdf)


def main():
    args = parse_args()
    result = build_composition(args.report_dir, write_pdf=not args.no_pdf)
    if result is None:
        raise SystemExit("No A/B/C candidates under {}; nothing to compose.".format(args.report_dir))
    html_path, _, pdf_path = result
    print("Wrote HTML:", html_path.resolve())
    if pdf_path is not None:
        print("Wrote PDF:", pdf_path.resolve())


if __name__ == "__main__":
    main()
