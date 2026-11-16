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
        span_left = entry.get("span_left", entry["x1"])
        span_right = entry.get("span_right", entry["x2"])
        placed = False
        for lane_index, last_end in enumerate(lane_last_end):
            if span_left > last_end + pad_px:
                entry["lane"] = lane_index
                lane_last_end[lane_index] = span_right
                placed = True
                break
        if not placed:
            entry["lane"] = len(lane_last_end)
            lane_last_end.append(span_right)
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


def _clamp_center_x(mid_x, label_text, font_size, min_x, max_x, bold=False):
    half_w = rc.text_width(str(label_text), font_size, bold=bold) / 2.0
    if max_x - min_x <= 2.0 * half_w:
        return (min_x + max_x) / 2.0
    return max(min_x + half_w, min(max_x - half_w, float(mid_x)))


def _reference_ticks(ref_len):
    if ref_len <= 1:
        return [1]
    ticks = [1]
    for frac in (0.25, 0.50, 0.75):
        val = int(round(1 + (ref_len - 1) * frac))
        if 1 < val < ref_len and val not in ticks:
            ticks.append(val)
    if ref_len not in ticks:
        ticks.append(ref_len)
    return ticks


def build_html(reference, reference_label, grouped_hits):
    """Build the Nature-style two-panel nucleotide alignment figure via rc.Figure."""
    reference_length = next(iter(next(iter(grouped_hits.values())))).slen if grouped_hits else 0
    left_margin = 420
    right_margin = 80
    track_width = 980
    fig_width = left_margin + track_width + right_margin
    label_col_x = 28
    label_col_max_w = left_margin - label_col_x - 24
    baseline_height = 12
    hit_height = baseline_height
    lane_step = 24

    title_style = {"size": 13.5, "fill": "#222222", "bold": True}
    panel_tag_style = {"size": 15, "fill": "#111111", "bold": True}
    panel_title_style = {"size": 12, "fill": "#222222", "bold": True}
    label_style = {"size": 10.5, "fill": "#222222", "bold": True}
    subtle_style = {"size": 9.5, "fill": "#555555"}
    axis_style = {"size": 9.5, "fill": "#444444"}
    legend_style = {"size": 10, "fill": "#444444"}
    coord_style = {"size": 9.5, "fill": "#333333", "anchor": "middle"}

    sorted_items = sorted(
        grouped_hits.items(),
        key=lambda item: (primary_reference_start(item[1]), -best_bitscore(item[1]), -max(hit.qlen for hit in item[1])),
    )

    max_query_length = max((max(hit.qlen for hit in hits) for hits in grouped_hits.values()), default=0)
    total_hsps = sum(len(hits) for _, hits in sorted_items)
    contig_positions = {}
    axis_min = 1
    axis_max = max(1, reference_length)
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
    track_right_limit = fig_width - 16
    reference_entries = []
    for hit in all_reference_hits:
        ref_hit_x1 = scaled_x(hit.subject_min, axis_min, axis_span, ref_x0, track_width)
        ref_hit_width = scaled_width(hit.subject_max - hit.subject_min + 1, axis_span, track_width)
        ref_hit_x2 = ref_hit_x1 + ref_hit_width
        s_label = "s:{}-{}".format(hit.subject_min, hit.subject_max)
        label_cx = _clamp_center_x(ref_hit_x1 + ref_hit_width / 2.0, s_label, 9.5, ref_x0, track_right_limit)
        half_lw = rc.text_width(s_label, 9.5) / 2.0
        reference_entries.append(
            {
                "hit": hit,
                "x1": ref_hit_x1,
                "x2": ref_hit_x2,
                "width": ref_hit_width,
                "label": s_label,
                "label_cx": label_cx,
                "span_left": min(ref_hit_x1, label_cx - half_lw) - 3.0,
                "span_right": max(ref_hit_x2, label_cx + half_lw) + 3.0,
            }
        )
    lane_count = assign_reference_label_lanes(reference_entries, 10.0)

    # Pre-layout header lines
    header_max_w = fig_width - 360
    title_lines = rc.wrap_text(reference_label or reference, header_max_w, size=13.5, bold=True)
    meta_text = (
        "Reference length: {} bp | Longest contig: {} bp | Displayed contigs: {} | Displayed HSPs: {} | "
        "Display axis: {}..{} bp | pident range: {:.3f}-{:.3f}"
    ).format(
        reference_length,
        max_query_length,
        len(sorted_items),
        total_hsps,
        axis_min,
        axis_max,
        identity_low,
        identity_high,
    )
    meta_lines = rc.wrap_text(meta_text, header_max_w, size=9.5)
    header_bottom = max(78, 28 + len(title_lines) * 17 + len(meta_lines) * 14 + 14)

    # Panel a geometry
    panel_a_title_y = header_bottom + 24
    first_lane_top_y = panel_a_title_y + 30
    ref_bar_top = first_lane_top_y + max(0, lane_count - 1) * lane_step + 14
    ref_y = ref_bar_top + baseline_height / 2.0
    ruler_y = ref_y + baseline_height / 2.0 + 8

    # Pre-layout Panel b rows so figure height is exact
    panel_b_title_y = ruler_y + 42
    rows_start_y = panel_b_title_y + 34
    row_layouts = []
    cursor_y = rows_start_y
    for qseqid, hits in sorted_items:
        qlen = contig_positions[qseqid]["qlen"]
        contig_start_abs = contig_positions[qseqid]["contig_start_abs"]
        contig_orientation = contig_positions[qseqid]["orientation"]
        best = max(hit.bitscore for hit in hits)
        anchor = contig_positions[qseqid]["anchor"]

        id_lines = rc.wrap_text(qseqid, label_col_max_w, size=10.5, bold=True)
        sub_line_1 = "len={} bp | HSPs={} | bitscore={:.1f} | strand={}".format(
            qlen, len(hits), best, contig_orientation
        )
        sub_line_2 = "anchor: q:{}-{} -> s:{}-{}".format(
            anchor.query_min, anchor.query_max, anchor.subject_min, anchor.subject_max
        )
        sub_lines = (
            rc.wrap_text(sub_line_1, label_col_max_w, size=9.5)
            + rc.wrap_text(sub_line_2, label_col_max_w, size=9.5)
        )

        q_entries = []
        for hit in sorted(hits, key=lambda h: (h.query_min, h.query_max, -h.bitscore)):
            qry_hit_start_abs, qry_hit_end_abs = query_interval_absolute(
                hit, contig_start_abs, qlen, contig_orientation
            )
            qry_hit_x1 = scaled_x(qry_hit_start_abs, axis_min, axis_span, ref_x0, track_width)
            qry_hit_width = scaled_width(qry_hit_end_abs - qry_hit_start_abs + 1, axis_span, track_width)
            q_label = "q:{}-{}".format(hit.query_min, hit.query_max)
            q_cx = _clamp_center_x(qry_hit_x1 + qry_hit_width / 2.0, q_label, 9.5, ref_x0, track_right_limit)
            half_qw = rc.text_width(q_label, 9.5) / 2.0
            q_entries.append(
                {
                    "hit": hit,
                    "x1": qry_hit_x1,
                    "width": qry_hit_width,
                    "label": q_label,
                    "label_cx": q_cx,
                    "span_left": q_cx - half_qw - 2.0,
                    "span_right": q_cx + half_qw + 2.0,
                }
            )
        q_lanes = rc.assign_text_lanes(q_entries, pad_px=10.0)
        left_block_h = len(id_lines) * 13 + len(sub_lines) * 12 + 8
        track_block_h = 34 + max(1, q_lanes) * 15
        row_h = max(62, left_block_h + 14, track_block_h + 12)
        row_layouts.append(
            {
                "qseqid": qseqid,
                "qlen": qlen,
                "contig_start_abs": contig_start_abs,
                "contig_orientation": contig_orientation,
                "id_lines": id_lines,
                "sub_lines": sub_lines,
                "q_entries": q_entries,
                "row_top": cursor_y,
                "row_h": row_h,
            }
        )
        cursor_y += row_h

    total_height = cursor_y + 28
    fig = rc.Figure(fig_width, total_height)

    # Draw Header
    ty = 26
    for line in title_lines:
        fig.text(label_col_x, ty, line, **title_style)
        ty += 17
    ty += 2
    for line in meta_lines:
        fig.text(label_col_x, ty, line, **subtle_style)
        ty += 14

    legend_w = 230
    legend_x = fig_width - right_margin - legend_w
    legend_y = 32
    legend_h = 11
    fig.text(
        legend_x + legend_w / 2.0,
        legend_y - 8,
        "Nucleotide identity (pident)",
        size=10,
        fill="#333333",
        anchor="middle",
        bold=True,
    )
    fig.gradient_rect(
        legend_x,
        legend_y,
        legend_w,
        legend_h,
        lambda frac: pident_to_color(identity_low + (identity_high - identity_low) * frac, identity_low, identity_high),
        stroke=None,
    )
    fig.text(legend_x, legend_y + legend_h + 14, "{:.1f}%".format(identity_low), anchor="start", **legend_style)
    fig.text(
        legend_x + legend_w / 2.0,
        legend_y + legend_h + 14,
        "{:.1f}%".format(0.5 * (identity_low + identity_high)),
        anchor="middle",
        **legend_style
    )
    fig.text(legend_x + legend_w, legend_y + legend_h + 14, "{:.1f}%".format(identity_high), anchor="end", **legend_style)

    fig.polyline([(label_col_x, header_bottom), (fig_width - 28, header_bottom)], stroke="#E2E6EA", stroke_width=1.0)

    # Draw Panel a: Reference track & non-overlapping HSP interval callouts
    fig.text(label_col_x, panel_a_title_y, "a", **panel_tag_style)
    fig.text(label_col_x + 20, panel_a_title_y, "Reference genome & aligned HSP intervals", **panel_title_style)

    ref_Lines = rc.wrap_text(reference, label_col_max_w, size=10.5, bold=True)
    for idx, rline in enumerate(ref_Lines):
        fig.text(label_col_x, ref_y - 1 + idx * 13, rline, **label_style)
    fig.text(
        label_col_x,
        ref_y - 1 + len(ref_Lines) * 13,
        "Reference scaffold ({} bp)".format(reference_length),
        **subtle_style
    )

    ref_start_x = scaled_x(1, axis_min, axis_span, ref_x0, track_width)
    ref_bar_width = scaled_width(reference_length, axis_span, track_width)
    ref_end_x = ref_start_x + ref_bar_width
    fig.rect(
        ref_start_x,
        ref_y - baseline_height / 2.0,
        ref_bar_width,
        baseline_height,
        rx=2,
        fill="#272727",
        stroke=None,
    )

    # Coordinate ruler under reference bar
    fig.polyline([(ref_start_x, ruler_y), (ref_end_x, ruler_y)], stroke="#555555", stroke_width=1.0)
    last_tick_right = -1e9
    for tick_pos in _reference_ticks(reference_length):
        tx = scaled_x(tick_pos, axis_min, axis_span, ref_x0, track_width)
        fig.polyline([(tx, ruler_y), (tx, ruler_y + 4.0)], stroke="#555555", stroke_width=1.0)
        t_label = "{} bp".format(tick_pos) if tick_pos in (1, reference_length) else str(tick_pos)
        t_cx = _clamp_center_x(tx, t_label, 9.5, ref_x0, track_right_limit)
        half_tw = rc.text_width(t_label, 9.5) / 2.0
        if t_cx - half_tw > last_tick_right + 10.0:
            fig.text(t_cx, ruler_y + 17.0, t_label, anchor="middle", **axis_style)
            last_tick_right = t_cx + half_tw

    # Self-contained lane interval brackets above the reference bar
    for entry in reference_entries:
        hit = entry["hit"]
        ref_hit_x1 = entry["x1"]
        ref_hit_x2 = entry["x2"]
        lane = entry["lane"]
        lane_marker_y = ref_bar_top - 10.0 - lane * lane_step
        fig.polyline(
            [
                (ref_hit_x1, lane_marker_y + 3.0),
                (ref_hit_x1, lane_marker_y),
                (ref_hit_x2, lane_marker_y),
                (ref_hit_x2, lane_marker_y + 3.0),
            ],
            stroke="#272727",
            stroke_width=1.2,
            title=title_text(hit),
        )
        fig.text(
            entry["label_cx"],
            lane_marker_y - 6.0,
            entry["label"],
            **coord_style
        )

    # Draw Panel b: Assembled contig tracks
    fig.text(label_col_x, panel_b_title_y, "b", **panel_tag_style)
    fig.text(
        label_col_x + 20,
        panel_b_title_y,
        "Assembled contig alignments (ordered by reference position & bitscore)",
        **panel_title_style
    )

    for layout in row_layouts:
        qlen = layout["qlen"]
        contig_start_abs = layout["contig_start_abs"]
        contig_orientation = layout["contig_orientation"]
        row_top = layout["row_top"]
        query_y = row_top + 18.0

        # Left column text (strictly inside label_col_max_w)
        ly = row_top + 10.0
        for line in layout["id_lines"]:
            fig.text(label_col_x, ly, line, **label_style)
            ly += 13.0
        for line in layout["sub_lines"]:
            fig.text(label_col_x, ly, line, **subtle_style)
            ly += 12.0

        # Contig scaffold bar
        query_start_x = scaled_x(contig_start_abs, axis_min, axis_span, ref_x0, track_width)
        query_bar_width = scaled_width(qlen, axis_span, track_width)
        query_end_x = query_start_x + query_bar_width
        fig.rect(
            query_start_x,
            query_y - baseline_height / 2.0,
            query_bar_width,
            baseline_height,
            rx=2,
            fill="#DCE1E7",
            stroke=None,
        )

        # Contig endpoint / orientation labels placed cleanly ABOVE the contig bar
        left_coord = qlen if contig_orientation == "-" else 1
        right_coord = 1 if contig_orientation == "-" else qlen
        left_str = "{} bp ({})".format(left_coord, contig_orientation)
        right_str = "{} bp".format(right_coord)
        left_w = rc.text_width(left_str, 9.5)
        right_w = rc.text_width(right_str, 9.5)
        lx_start = max(ref_x0, min(track_right_limit - left_w, query_start_x))
        rx_end = max(ref_x0 + right_w, min(track_right_limit, query_end_x))
        if query_bar_width >= left_w + right_w + 18.0 and (rx_end - right_w) > (lx_start + left_w + 10.0):
            fig.text(lx_start, query_y - 10.0, left_str, anchor="start", **axis_style)
            fig.text(rx_end, query_y - 10.0, right_str, anchor="end", **axis_style)
        else:
            combined = "{}..{} bp ({})".format(left_coord, right_coord, contig_orientation)
            c_cx = _clamp_center_x(query_start_x + query_bar_width / 2.0, combined, 9.5, ref_x0, track_right_limit)
            fig.text(c_cx, query_y - 10.0, combined, anchor="middle", **axis_style)

        # Coloured HSP blocks & lane-staggered q:start-end labels BELOW the contig bar
        for q_entry in layout["q_entries"]:
            hit = q_entry["hit"]
            color = pident_to_color(hit.pident, identity_low, identity_high)
            fig.rect(
                q_entry["x1"],
                query_y - hit_height / 2.0,
                q_entry["width"],
                hit_height,
                rx=2,
                fill=color,
                stroke=None,
                title=title_text(hit),
            )
            label_y = query_y + baseline_height / 2.0 + 13.0 + q_entry["lane"] * 14.0
            fig.text(
                q_entry["label_cx"],
                label_y,
                q_entry["label"],
                **coord_style
            )

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
