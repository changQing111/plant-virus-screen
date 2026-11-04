#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path

_LEGACY_DIR = str(Path(__file__).resolve().parent)
if _LEGACY_DIR not in sys.path:
    sys.path.insert(0, _LEGACY_DIR)

import render_common as rc


class BlastnHit(object):
    __slots__ = (
        "qseqid",
        "sseqid",
        "qlen",
        "slen",
        "pident",
        "length",
        "qcovs",
        "qcovhsp",
        "qcovus",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    )

    def __init__(self, row):
        self.qseqid = row[0]
        self.sseqid = row[1]
        self.qlen = int(row[2])
        self.slen = int(row[3])
        self.pident = float(row[4])
        self.length = int(row[5])
        self.qcovs = float(row[6])
        self.qcovhsp = float(row[7])
        self.qcovus = float(row[8])
        self.mismatch = int(row[9])
        self.gapopen = int(row[10])
        self.qstart = int(row[11])
        self.qend = int(row[12])
        self.sstart = int(row[13])
        self.send = int(row[14])
        self.evalue = float(row[15])
        self.bitscore = float(row[16])

    @property
    def query_min(self):
        return min(self.qstart, self.qend)

    @property
    def query_max(self):
        return max(self.qstart, self.qend)

    @property
    def subject_min(self):
        return min(self.sstart, self.send)

    @property
    def subject_max(self):
        return max(self.sstart, self.send)

    @property
    def strand(self):
        if (self.qend - self.qstart) * (self.send - self.sstart) >= 0:
            return "+"
        return "-"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize blastn hits for one reference accession as a self-contained HTML summary with a native vector PDF export."
    )
    parser.add_argument("--blastn", required=True, help="blastn result with the 17-column format used by run_blastn.sh")
    parser.add_argument("--reference", required=True, help="Reference accession to visualize, e.g. NC_030378.1")
    parser.add_argument("--out", required=True, help="Output HTML file")
    parser.add_argument(
        "--final-tsv",
        help="Optional final_retained_candidates.tsv. If set, only contigs whose best_nt_hit starts with --reference are plotted.",
    )
    parser.add_argument(
        "--top-contigs",
        type=int,
        default=0,
        help="Optional limit on the number of contigs to display after sorting by best bitscore (0 means all).",
    )
    parser.add_argument(
        "--min-align-len",
        type=int,
        default=0,
        help="Optional minimum HSP alignment length to display.",
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


def load_blastn_hits(path, reference, min_align_len):
    hits = []
    with Path(path).open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[1] != reference:
                continue
            hit = BlastnHit(row)
            if hit.length < min_align_len:
                continue
            hits.append(hit)
    return hits


def load_allowed_contigs(final_tsv, reference):
    allowed = set()
    with Path(final_tsv).open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            best_nt_hit = (row.get("best_nt_hit") or "").strip()
            if best_nt_hit.startswith(reference):
                allowed.add((row.get("representative_query_id") or "").strip())
    return allowed


def load_reference_label(final_tsv, reference):
    if not final_tsv:
        return reference
    with Path(final_tsv).open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            best_nt_hit = (row.get("best_nt_hit") or "").strip()
            if best_nt_hit.startswith(reference):
                return best_nt_hit
    return reference


def group_hits_by_contig(hits):
    grouped = {}
    for hit in hits:
        grouped.setdefault(hit.qseqid, []).append(hit)
    return grouped


def best_bitscore(hit_list):
    return max(hit.bitscore for hit in hit_list)


def anchor_hit(hit_list):
    return max(hit_list, key=lambda hit: (hit.bitscore, hit.length, hit.pident))


def primary_reference_start(hit_list):
    return anchor_hit(hit_list).subject_min


def observed_identity_range(grouped_hits):
    return rc.observed_identity_range(grouped_hits)


def pident_to_color(pident, low, high):
    return rc.pident_to_color(pident, low, high)


def scaled_x(position, axis_min, axis_span, x0, width):
    if axis_span <= 0:
        return x0
    return x0 + (float(position - axis_min) / float(axis_span)) * width


def scaled_width(span_length, axis_span, width):
    if axis_span <= 0:
        return 0.0
    return max(1.5, (float(span_length) / float(axis_span)) * width)


def title_text(hit):
    return (
        "contig={qseqid} | q={qstart}-{qend}/{qlen} | s={sstart}-{send}/{slen} | "
        "len={length} | pident={pident:.3f} | qcovs={qcovs:.1f} | bitscore={bitscore:.1f} | strand={strand}"
    ).format(
        qseqid=hit.qseqid,
        qstart=hit.qstart,
        qend=hit.qend,
        qlen=hit.qlen,
        sstart=hit.sstart,
        send=hit.send,
        slen=hit.slen,
        length=hit.length,
        pident=hit.pident,
        qcovs=hit.qcovs,
        bitscore=hit.bitscore,
        strand=hit.strand,
    )


def contig_absolute_bounds(anchor, qlen):
    if anchor.strand == "+":
        start = anchor.subject_min - anchor.query_min + 1
    else:
        start = anchor.subject_min + anchor.query_max - qlen
    end = start + qlen - 1
    return start, end


def query_interval_absolute(hit, contig_start_abs, qlen, contig_orientation):
    if contig_orientation == "+":
        start = contig_start_abs + hit.query_min - 1
        end = contig_start_abs + hit.query_max - 1
    else:
        start = contig_start_abs + (qlen - hit.query_max)
        end = contig_start_abs + (qlen - hit.query_min)
    return start, end


def assign_reference_label_lanes(entries, pad_px):
    lane_last_end = []
    for entry in entries:
        placed = False
        for lane_index, last_end in enumerate(lane_last_end):
            if entry["x1"] > last_end + pad_px:
                entry["lane"] = lane_index
                lane_last_end[lane_index] = entry["x2"]
                placed = True
                break
        if not placed:
            entry["lane"] = len(lane_last_end)
            lane_last_end.append(entry["x2"])
    return len(lane_last_end)


def deduplicate_reference_hits(hits):
    seen = set()
    unique = []
    for hit in hits:
        key = (hit.subject_min, hit.subject_max)
        if key in seen:
            continue
        seen.add(key)
        unique.append(hit)
    return unique


def build_html(reference, reference_label, grouped_hits):
    """Build the figure via the shared rendering layer.

    Layout and visual conventions are unchanged from the previous
    per-script SVG implementation; drawing goes through rc.Figure so the
    same figure is rendered to inline SVG (HTML) and native vector PDF.
    """
    reference_length = next(iter(next(iter(grouped_hits.values())))).slen if grouped_hits else 0
    left_margin = 220
    right_margin = 80
    track_width = 1100
    row_height = 54
    baseline_height = 14
    hit_height = baseline_height
    gap_after_reference = 42
    lane_step = 14

    label_style = {"size": 12, "fill": "#222222"}
    subtle_style = {"size": 11, "fill": "#666666"}
    axis_style = {"size": 11, "fill": "#555555"}
    legend_style = {"size": 11, "fill": "#444444"}
    coord_style = {"size": 9, "fill": "#333333", "anchor": "middle"}

    sorted_items = sorted(
        grouped_hits.items(),
        key=lambda item: (primary_reference_start(item[1]), -best_bitscore(item[1]), -max(hit.qlen for hit in item[1])),
    )

    max_query_length = max((max(hit.qlen for hit in hits) for hits in grouped_hits.values()), default=0)
    contig_positions = {}
    axis_min = 1
    axis_max = reference_length
    for qseqid, hits in sorted_items:
        qlen = max(hit.qlen for hit in hits)
        anchor = anchor_hit(hits)
        contig_start_abs, contig_end_abs = contig_absolute_bounds(anchor, qlen)
        contig_positions[qseqid] = {
            "anchor": anchor,
            "qlen": qlen,
            "contig_start_abs": contig_start_abs,
            "contig_end_abs": contig_end_abs,
            "orientation": anchor.strand,
        }
        axis_min = min(axis_min, contig_start_abs)
        axis_max = max(axis_max, contig_end_abs)
    axis_span = max(1, axis_max - axis_min + 1)

    identity_low, identity_high = observed_identity_range(grouped_hits)

    all_reference_hits = []
    for hits in grouped_hits.values():
        all_reference_hits.extend(hits)
    all_reference_hits = deduplicate_reference_hits(all_reference_hits)
    all_reference_hits.sort(key=lambda hit: (hit.subject_min, hit.subject_max, -hit.bitscore))

    ref_x0 = left_margin
    reference_entries = []
    for hit in all_reference_hits:
        ref_hit_x1 = scaled_x(hit.subject_min, axis_min, axis_span, ref_x0, track_width)
        ref_hit_width = scaled_width(hit.subject_max - hit.subject_min + 1, axis_span, track_width)
        reference_entries.append(
            {
                "hit": hit,
                "x1": ref_hit_x1,
                "x2": ref_hit_x1 + ref_hit_width,
                "width": ref_hit_width,
            }
        )
    lane_count = assign_reference_label_lanes(reference_entries, 10.0)

    ref_y = 90 + max(0, lane_count - 1) * lane_step
    header_height = ref_y + gap_after_reference
    total_height = header_height + len(sorted_items) * row_height + 40

    fig = rc.Figure(left_margin + track_width + right_margin, total_height)

    meta_text = (
        "Reference length: {} bp | Longest contig: {} bp | Displayed contigs: {} | Displayed HSPs: {} | "
        "Display axis: {}..{} bp | "
        "pident range: {:.3f}-{:.3f}"
    ).format(
        reference_length,
        max_query_length,
        len(sorted_items),
        sum(len(hits) for _, hits in sorted_items),
        axis_min,
        axis_max,
        identity_low,
        identity_high,
    )

    fig.text(20, 28, "Reference track in black, contig tracks in gray. Colored blocks mark HSPs.", **subtle_style)

    legend_x = left_margin + track_width - 260
    legend_y = 30
    legend_w = 220
    legend_h = 12
    fig.gradient_rect(
        legend_x,
        legend_y,
        legend_w,
        legend_h,
        lambda frac: pident_to_color(identity_low + (identity_high - identity_low) * frac, identity_low, identity_high),
        stroke="#bbbbbb",
    )
    fig.text(legend_x, legend_y + 28, "{:.2f}%".format(identity_low), **legend_style)
    fig.text(legend_x + legend_w - 36, legend_y + 28, "{:.2f}%".format(identity_high), **legend_style)
    fig.text(legend_x + legend_w / 2 - 18, legend_y + 28, "pident", **legend_style)

    ref_marker_y = ref_y - baseline_height / 2 - 7
    fig.text(20, ref_y + 4, reference, **label_style)
    ref_start_x = scaled_x(1, axis_min, axis_span, ref_x0, track_width)
    ref_bar_width = scaled_width(reference_length, axis_span, track_width)
    fig.rect(
        ref_start_x,
        ref_y - baseline_height / 2,
        ref_bar_width,
        baseline_height,
        rx=2,
        fill="#111111",
        stroke="#111111",
        stroke_width=1.0,
    )
    fig.text(ref_start_x - 8, ref_y + 24, "1", **axis_style)
    ref_end_x = ref_start_x + ref_bar_width
    fig.text(ref_end_x - 30, ref_y + 24, reference_length, **axis_style)

    # Keep the reference bar pure black, but mark each HSP's true interval with
    # a black bracket above the bar so the subject coordinates remain visible.
    for entry in reference_entries:
        hit = entry["hit"]
        ref_hit_x1 = entry["x1"]
        ref_hit_x2 = entry["x2"]
        ref_mid_x = ref_hit_x1 + entry["width"] / 2.0
        lane = entry["lane"]
        lane_marker_y = ref_marker_y - lane * lane_step
        fig.polyline(
            [
                (ref_hit_x1, lane_marker_y),
                (ref_hit_x1, ref_y - baseline_height / 2),
                (ref_hit_x2, ref_y - baseline_height / 2),
                (ref_hit_x2, lane_marker_y),
            ],
            stroke="#111111",
            stroke_width=1.5,
            title=title_text(hit),
        )
        fig.text(
            ref_mid_x,
            lane_marker_y - 3,
            "s:{}-{}".format(hit.subject_min, hit.subject_max),
            **coord_style
        )

    current_y = ref_y + gap_after_reference
    for qseqid, hits in sorted_items:
        qlen = contig_positions[qseqid]["qlen"]
        contig_start_abs = contig_positions[qseqid]["contig_start_abs"]
        contig_orientation = contig_positions[qseqid]["orientation"]
        best = max(hit.bitscore for hit in hits)
        label_main = "{}".format(qseqid)
        anchor = contig_positions[qseqid]["anchor"]
        label_sub = "len={} | HSPs={} | best bitscore={:.1f} | strand={} | anchored by q:{}-{} -> s:{}-{}".format(
            qlen,
            len(hits),
            best,
            contig_orientation,
            anchor.query_min,
            anchor.query_max,
            anchor.subject_min,
            anchor.subject_max,
        )
        fig.text(20, current_y + 4, label_main, **label_style)
        fig.text(20, current_y + 18, label_sub, **subtle_style)

        query_y = current_y
        query_start_x = scaled_x(contig_start_abs, axis_min, axis_span, ref_x0, track_width)
        query_bar_width = scaled_width(qlen, axis_span, track_width)
        fig.rect(
            query_start_x,
            query_y - baseline_height / 2,
            query_bar_width,
            baseline_height,
            rx=2,
            fill="#a0a0a0",
            stroke="#8d8d8d",
            stroke_width=1.0,
        )
        left_label = qlen if contig_orientation == "-" else 1
        right_label = 1 if contig_orientation == "-" else qlen
        fig.text(query_start_x - 16, query_y + 24, left_label, **axis_style)
        query_end_x = query_start_x + query_bar_width
        fig.text(query_end_x - 30, query_y + 24, right_label, **axis_style)

        for hit in hits:
            color = pident_to_color(hit.pident, identity_low, identity_high)
            qry_hit_start_abs, qry_hit_end_abs = query_interval_absolute(
                hit, contig_start_abs, qlen, contig_orientation
            )
            qry_hit_x1 = scaled_x(qry_hit_start_abs, axis_min, axis_span, ref_x0, track_width)
            qry_hit_width = scaled_width(qry_hit_end_abs - qry_hit_start_abs + 1, axis_span, track_width)
            qry_mid_x = qry_hit_x1 + qry_hit_width / 2.0

            fig.rect(
                qry_hit_x1,
                query_y - hit_height / 2,
                qry_hit_width,
                hit_height,
                rx=2,
                fill=color,
                title=title_text(hit),
            )
            fig.text(
                qry_mid_x,
                query_y + baseline_height / 2 + 12,
                "q:{}-{}".format(hit.query_min, hit.query_max),
                **coord_style
            )

        current_y += row_height

    return rc.render_html(fig, reference_label + " blastn summary", reference_label, meta_text)


def write_output_bundle(out_path, html_doc, write_svg=True, write_pdf=True):
    """Delegate to the shared rendering layer (see render_common)."""
    return rc.write_output_bundle(out_path, html_doc, write_svg=write_svg, write_pdf=write_pdf)


def main():
    args = parse_args()
    hits = load_blastn_hits(args.blastn, args.reference, args.min_align_len)
    if not hits:
        raise SystemExit("No blastn hits found for reference {}".format(args.reference))

    if args.final_tsv:
        allowed = load_allowed_contigs(args.final_tsv, args.reference)
        hits = [hit for hit in hits if hit.qseqid in allowed]
        if not hits:
            raise SystemExit("No blastn hits remain after filtering by final candidates for {}".format(args.reference))

    grouped = group_hits_by_contig(hits)
    sorted_names = sorted(grouped, key=lambda name: (best_bitscore(grouped[name]), max(hit.qlen for hit in grouped[name])), reverse=True)
    if args.top_contigs > 0:
        sorted_names = sorted_names[: args.top_contigs]
    grouped = {name: grouped[name] for name in sorted_names}

    reference_label = load_reference_label(args.final_tsv, args.reference)
    html_doc = build_html(args.reference, reference_label, grouped)
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
