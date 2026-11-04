#!/usr/bin/env python3
"""Coverage depth curve for one candidate contig.

Draws the per-position read depth (samtools depth -a -Q 20 output) along
the contig as a binned curve with an area fill, a dashed mean-depth
reference line, and MAPQ>=20 breadth statistics. Built on the shared
rendering layer, so each figure is written as self-contained HTML (inline
SVG, hover tooltips) plus a native vector PDF.

Python 3.6+ compatible, standard library only.
"""

import argparse
import csv
import sys
from pathlib import Path

_LEGACY_DIR = str(Path(__file__).resolve().parent)
if _LEGACY_DIR not in sys.path:
    sys.path.insert(0, _LEGACY_DIR)

import render_common as rc


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize MAPQ>=20 per-position depth for one candidate contig as HTML + native vector PDF."
    )
    parser.add_argument("--depth", required=True, help="samtools depth -a output (contig\\tposition\\tdepth), MAPQ-filtered")
    parser.add_argument("--contig", required=True, help="Candidate contig id to visualize")
    parser.add_argument("--out", required=True, help="Output HTML file")
    parser.add_argument(
        "--final-tsv",
        help="Optional plot/final candidates TSV used to annotate grade and best hits.",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=1000,
        help="Maximum number of curve bins along the contig (default 1000).",
    )
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="Only write HTML and skip automatic PDF export.",
    )
    parser.add_argument(
        "--no-svg",
        action="store_true",
        help="Deprecated no-op kept for compatibility: standalone SVG output was replaced by native vector PDF.",
    )
    return parser.parse_args()


def load_depth_positions(path, contigs):
    """Single streaming scan of a samtools depth file.

    Returns {contig_id: [depth_at_position_1, depth_at_position_2, ...]}.
    Only the requested contigs are collected; missing positions stay 0
    only if samtools emitted them (``-a`` guarantees every position).
    """
    wanted = set(contigs)
    depths = {contig: [] for contig in wanted}
    with Path(path).open() as handle:
        for line in handle:
            parts = line.split("\t", 2)
            if len(parts) < 3:
                continue
            contig = parts[0]
            if contig not in wanted:
                continue
            try:
                position = int(parts[1])
                depth = int(parts[2])
            except ValueError:
                continue
            series = depths[contig]
            # Fill any gap (defensive: -a should emit every position).
            gap = position - 1 - len(series)
            if gap > 0:
                series.extend([0] * gap)
            series.append(depth)
    return depths


def bin_depths(values, bins):
    """Mean depth per bin; bin count never exceeds len(values) or bins."""
    count = len(values)
    if count == 0:
        return []
    bins = max(1, min(int(bins), count))
    bin_size = -(-count // bins)
    binned = []
    for start in range(0, count, bin_size):
        chunk = values[start : start + bin_size]
        binned.append(sum(chunk) / float(len(chunk)))
    return binned


def compute_stats(values):
    """Coverage statistics from the per-position depth series."""
    count = len(values)
    if count == 0:
        return {"length": 0, "mean": 0.0, "max": 0, "p99": 0.0, "breadth_1x": 0.0, "breadth_5x": 0.0}
    total = sum(values)
    ge1 = 0
    ge5 = 0
    for v in values:
        if v >= 1:
            ge1 += 1
        if v >= 5:
            ge5 += 1
    ordered = sorted(values)
    p99_index = min(count - 1, max(0, int(round(0.99 * (count - 1)))))
    return {
        "length": count,
        "mean": total / float(count),
        "max": ordered[-1],
        "p99": float(ordered[p99_index]),
        "breadth_1x": ge1 / float(count),
        "breadth_5x": ge5 / float(count),
    }


def load_annotation(final_tsv, contig):
    """Optional annotation row (grade, best hits) from the candidates TSV."""
    annotation = {
        "evidence_grade": "",
        "putative_virus": "",
        "blastn_anno": "",
        "blastx_anno": "",
        "mean_depth_tsv": "",
        "breadth_1x_tsv": "",
    }
    if not final_tsv:
        return annotation
    with Path(final_tsv).open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row.get("query_id") or row.get("contig_id") or "") == contig:
                annotation["evidence_grade"] = (row.get("evidence_grade") or "").strip()
                annotation["putative_virus"] = (
                    row.get("putative_virus") or row.get("candidate_name") or ""
                ).strip()
                annotation["blastn_anno"] = (row.get("blastn_anno") or "").strip()
                annotation["blastx_anno"] = (row.get("blastx_anno") or "").strip()
                annotation["mean_depth_tsv"] = (row.get("mapq20_mean_depth") or "").strip()
                annotation["breadth_1x_tsv"] = (row.get("mapq20_breadth_1x") or "").strip()
                break
    return annotation


def build_html(contig_id, depth_values, annotation, bins=1000):
    """Coverage curve figure via the shared rendering layer."""
    label_style = {"size": 12, "fill": "#222222"}
    subtle_style = {"size": 11, "fill": "#666666"}
    axis_style = {"size": 11, "fill": "#555555"}
    small_style = {"size": 10, "fill": "#555555"}

    stats = compute_stats(depth_values)
    binned = bin_depths(depth_values, bins)

    width = 1400
    x0 = 90.0
    x1 = 1340.0
    plot_y0 = 100.0
    plot_y1 = 340.0
    plot_h = plot_y1 - plot_y0
    total_height = 420

    # Robust y-axis cap: P99-based so a single spike cannot flatten the curve.
    y_cap = max(stats["p99"] * 1.25, stats["mean"] * 2.0, 1.0)
    capped = stats["max"] > y_cap

    fig = rc.Figure(width, total_height)

    header_bits = []
    if annotation.get("putative_virus"):
        header_bits.append(annotation["putative_virus"])
    if annotation.get("evidence_grade"):
        header_bits.append("grade " + annotation["evidence_grade"])
    fig.text(24, 34, contig_id, bold=True, **{"size": 13, "fill": "#222222"})
    if header_bits:
        fig.text(24, 54, " | ".join(header_bits), **subtle_style)

    meta_bits = [
        "Contig length: {} bp".format(stats["length"]),
        "mean depth (MAPQ>=20): {:.2f}".format(stats["mean"]),
        "breadth >=1x: {:.4f}".format(stats["breadth_1x"]),
        "breadth >=5x: {:.4f}".format(stats["breadth_5x"]),
        "true max depth: {}".format(stats["max"]),
    ]
    if capped:
        meta_bits.append("y-axis capped at {:.1f} (P99 x 1.25)".format(y_cap))
    meta_text = " | ".join(meta_bits)

    # Plot frame and baseline.
    fig.rect(x0, plot_y0, x1 - x0, plot_y1 - plot_y0, fill="#fcfcfc", stroke="#d0d0d0", stroke_width=1.0)

    def y_of(depth):
        clipped = min(float(depth), y_cap)
        return plot_y1 - (clipped / y_cap) * plot_h

    def x_of(frac):
        return x0 + frac * (x1 - x0)

    # Area fill + curve from binned depths.
    if binned:
        count = len(binned)
        step = 1.0 / count
        curve = [(x_of((index + 0.5) * step), y_of(value)) for index, value in enumerate(binned)]
        area = [(x0, plot_y1)] + curve + [(x1, plot_y1)]
        fig.polygon(area, fill="#b5d4f4")
        fig.polyline(curve, stroke="#185fa5", stroke_width=1.5, title="binned MAPQ>=20 depth")

    # Dashed mean-depth reference line.
    mean_y = y_of(stats["mean"])
    fig.polyline([(x0, mean_y), (x1, mean_y)], stroke="#d85a30", stroke_width=1.0, dash=(6, 4))
    fig.text(x1 + 6, mean_y + 4, "mean {:.1f}".format(stats["mean"]), **small_style)

    # Axis labels and ticks.
    fig.text(x0 - 58, plot_y0 + 8, "depth (MAPQ>=20)", **axis_style)
    fig.text(x0 - 34, plot_y1 + 4, "0", **axis_style)
    fig.text(x0 - 34, plot_y0 + 4, "{:.1f}".format(y_cap), **axis_style)

    mid = stats["length"] // 2
    fig.text(x_of(0.0), plot_y1 + 20, "0", anchor="middle", **axis_style)
    fig.text(x_of(0.5), plot_y1 + 20, mid, anchor="middle", **axis_style)
    fig.text(x_of(1.0), plot_y1 + 20, stats["length"], anchor="middle", **axis_style)
    fig.text(x_of(0.5), plot_y1 + 38, "position (bp)", anchor="middle", **axis_style)

    # Best-hit annotations (short) under the plot.
    line_y = plot_y1 + 66
    if annotation.get("blastn_anno"):
        fig.text(24, line_y, "nt: " + annotation["blastn_anno"], **subtle_style)
        line_y += 16
    if annotation.get("blastx_anno"):
        fig.text(24, line_y, "protein: " + annotation["blastx_anno"], **subtle_style)
        line_y += 16
    if capped:
        fig.text(
            24,
            line_y,
            "Note: y-axis capped at P99 x 1.25; true max depth {} is off-scale.".format(stats["max"]),
            **subtle_style
        )

    fig.text(24, 76, "Read depth along the contig; dashed line marks the mean depth.", **subtle_style)
    return rc.render_html(fig, contig_id + " coverage summary", contig_id, meta_text)


def write_output_bundle(out_path, html_doc, write_svg=True, write_pdf=True):
    """Delegate to the shared rendering layer (see render_common)."""
    return rc.write_output_bundle(out_path, html_doc, write_svg=write_svg, write_pdf=write_pdf)


def main():
    args = parse_args()
    depths = load_depth_positions(args.depth, [args.contig])
    values = depths.get(args.contig) or []
    if not values:
        raise SystemExit("No depth rows found for contig {}".format(args.contig))
    annotation = load_annotation(args.final_tsv, args.contig)
    html_doc = build_html(args.contig, values, annotation, bins=args.bins)
    out_path, svg_path, pdf_path = write_output_bundle(
        args.out,
        html_doc,
        write_svg=not args.no_svg,
        write_pdf=not args.no_pdf,
    )
    print("Wrote HTML:", out_path.resolve())
    if svg_path is not None:
        print("Wrote SVG:", svg_path.resolve())
    if pdf_path is not None:
        print("Wrote PDF:", pdf_path.resolve())


if __name__ == "__main__":
    main()
