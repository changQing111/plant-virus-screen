#!/usr/bin/env python3

import argparse
import csv
import re
import sys
from pathlib import Path

_LEGACY_DIR = str(Path(__file__).resolve().parent)
if _LEGACY_DIR not in sys.path:
    sys.path.insert(0, _LEGACY_DIR)

import render_common as rc


class BlastxHit(object):
    __slots__ = (
        "qseqid",
        "sseqid",
        "qlen",
        "slen",
        "pident",
        "length",
        "qcovhsp",
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
        self.qcovhsp = float(row[6])
        self.mismatch = int(row[7])
        self.gapopen = int(row[8])
        self.qstart = int(row[9])
        self.qend = int(row[10])
        self.sstart = int(row[11])
        self.send = int(row[12])
        self.evalue = float(row[13])
        self.bitscore = float(row[14])

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
        return "+" if self.qend >= self.qstart else "-"

    @property
    def frame_label(self):
        if self.strand == "+":
            frame = ((self.qstart - 1) % 3) + 1
            return "+{}".format(frame)
        frame = ((self.qend - 1) % 3) + 1
        return "-{}".format(frame)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize blastx hits for one viral protein accession as a self-contained HTML summary with a native vector PDF export."
    )
    parser.add_argument("--blastx", required=True, help="diamond blastx result file")
    parser.add_argument("--protein", required=True, help="Protein accession to visualize, e.g. tr|A0A165RBH9|A0A165RBH9_9VIRU")
    parser.add_argument("--out", required=True, help="Output HTML path")
    parser.add_argument(
        "--final-tsv",
        help="Optional final_retained_candidates.tsv. If set, only contigs whose best_protein_hit starts with --protein are plotted.",
    )
    parser.add_argument(
        "--protein-fasta",
        default="/home/changq/my_data/Document/viral/uniprotkb_taxonomy_id_10239_2026_04_03.fasta",
        help="Protein FASTA used to resolve the protein accession label.",
    )
    parser.add_argument(
        "--top-contigs",
        type=int,
        default=0,
        help="Optional cap on the number of contigs displayed after sorting.",
    )
    parser.add_argument(
        "--min-align-aa",
        type=int,
        default=0,
        help="Optional minimum amino-acid alignment length to display.",
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


def sanitize_filename(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def load_blastx_hits(path, protein, min_align_aa):
    hits = []
    with Path(path).open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[1] != protein:
                continue
            hit = BlastxHit(row)
            if hit.length < min_align_aa:
                continue
            hits.append(hit)
    return hits


def load_allowed_contigs(final_tsv, protein):
    if not final_tsv:
        return None
    allowed = set()
    with Path(final_tsv).open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            best_protein_hit = (row.get("best_protein_hit") or "").strip()
            if best_protein_hit.startswith(protein):
                rep = (row.get("representative_query_id") or "").strip()
                if rep:
                    allowed.add(rep)
    return allowed


def load_protein_label(protein, protein_fasta, final_tsv):
    if final_tsv:
        with Path(final_tsv).open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                best_protein_hit = (row.get("best_protein_hit") or "").strip()
                if best_protein_hit.startswith(protein):
                    return best_protein_hit
    fasta_path = Path(protein_fasta)
    if fasta_path.exists():
        with fasta_path.open() as handle:
            for line in handle:
                if line.startswith(">") and protein in line:
                    return line[1:].strip()
    return protein


def group_hits_by_contig(hits):
    grouped = {}
    for hit in hits:
        grouped.setdefault(hit.qseqid, []).append(hit)
    return grouped


def best_bitscore(hit_list):
    return max(hit.bitscore for hit in hit_list)


def anchor_hit(hit_list):
    return max(hit_list, key=lambda hit: (hit.bitscore, hit.length, hit.pident))


def observed_identity_range(grouped_hits):
    return rc.observed_identity_range(grouped_hits)


def pident_to_color(pident, low, high):
    return rc.pident_to_color(pident, low, high)


def scale(position, length, x0, width):
    if length <= 0:
        return x0
    return x0 + (float(position - 1) / float(length)) * width


def scaled_width(span, length, width):
    if length <= 0:
        return 0.0
    return max(1.5, (float(span) / float(length)) * width)


def hit_title(hit):
    return (
        "contig={qseqid} | q={qstart}-{qend}/{qlen} nt | s={sstart}-{send}/{slen} aa | "
        "len={length} aa | pident={pident:.3f} | qcovhsp={qcovhsp:.1f} | bitscore={bitscore:.1f} | "
        "strand={strand} | frame={frame}"
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
        qcovhsp=hit.qcovhsp,
        bitscore=hit.bitscore,
        strand=hit.strand,
        frame=hit.frame_label,
    )


def _clamp_center_x(mid_x, label_text, font_size, min_x, max_x, bold=False):
    half_w = rc.text_width(str(label_text), font_size, bold=bold) / 2.0
    if max_x - min_x <= 2.0 * half_w:
        return (min_x + max_x) / 2.0
    return max(min_x + half_w, min(max_x - half_w, float(mid_x)))


def _axis_ticks(max_len):
    if max_len <= 1:
        return [1]
    ticks = [1]
    for frac in (0.25, 0.50, 0.75):
        val = int(round(1 + (max_len - 1) * frac))
        if 1 < val < max_len and val not in ticks:
            ticks.append(val)
    if max_len not in ticks:
        ticks.append(max_len)
    return ticks


def build_html(protein, protein_label, grouped_hits):
    """Build the Nature-style two-panel figure (contig nt / protein aa) via rc.Figure."""
    protein_length = next(iter(next(iter(grouped_hits.values())))).slen if grouped_hits else 0
    max_query_length = max((max(hit.qlen for hit in hits) for hits in grouped_hits.values()), default=0)

    sorted_items = sorted(
        grouped_hits.items(),
        key=lambda item: (anchor_hit(item[1]).subject_min, -best_bitscore(item[1]), -max(hit.qlen for hit in item[1])),
    )
    total_hsps = sum(len(hits) for _, hits in sorted_items)
    identity_low, identity_high = observed_identity_range(grouped_hits)

    title_style = {"size": 13.5, "fill": "#222222", "bold": True}
    panel_tag_style = {"size": 15, "fill": "#111111", "bold": True}
    panel_title_style = {"size": 12, "fill": "#222222", "bold": True}
    label_style = {"size": 10.5, "fill": "#222222", "bold": True}
    subtle_style = {"size": 9.5, "fill": "#555555"}
    axis_style = {"size": 9.5, "fill": "#444444"}
    legend_style = {"size": 10, "fill": "#444444"}
    coord_style = {"size": 9.5, "fill": "#333333", "anchor": "middle"}

    fig_width = 1520
    label_col_x = 28
    contig_x0 = 400
    label_col_max_w = contig_x0 - label_col_x - 28
    contig_width = 560
    contig_right_limit = contig_x0 + contig_width + 20
    protein_x0 = 1030
    protein_width = 430
    protein_right_limit = fig_width - 16
    baseline_height = 12
    hit_height = 12

    # Pre-layout header
    header_max_w = fig_width - 360
    title_lines = rc.wrap_text(protein_label or protein, header_max_w, size=13.5, bold=True)
    meta_text = (
        "Protein length: {} aa | Longest contig: {} nt | Displayed contigs: {} | Displayed HSPs: {} | "
        "pident range: {:.3f}-{:.3f}"
    ).format(
        protein_length,
        max_query_length,
        len(sorted_items),
        total_hsps,
        identity_low,
        identity_high,
    )
    meta_lines = rc.wrap_text(meta_text, header_max_w, size=9.5)
    header_bottom = max(78, 28 + len(title_lines) * 17 + len(meta_lines) * 14 + 14)

    panel_header_y = header_bottom + 24
    protein_track_y = panel_header_y + 24
    ruler_y = protein_track_y + baseline_height / 2.0 + 8
    rows_start_y = ruler_y + 38

    row_layouts = []
    cursor_y = rows_start_y
    for qseqid, hits in sorted_items:
        qlen = max(hit.qlen for hit in hits)
        best = best_bitscore(hits)
        anchor = anchor_hit(hits)

        id_lines = rc.wrap_text(qseqid, label_col_max_w, size=10.5, bold=True)
        sub_str = "len={} nt | HSPs={} | bitscore={:.1f} | strand={} | frame={}".format(
            qlen, len(hits), best, anchor.strand, anchor.frame_label
        )
        sub_lines = rc.wrap_text(sub_str, label_col_max_w, size=9.5)

        q_entries = []
        s_entries = []
        for hit in sorted(hits, key=lambda h: (h.query_min, h.query_max, -h.bitscore)):
            if anchor.strand == "+":
                q_block_start = hit.query_min
                q_block_end = hit.query_max
            else:
                q_block_start = qlen - hit.query_max + 1
                q_block_end = qlen - hit.query_min + 1

            qx = scale(q_block_start, max_query_length, contig_x0, contig_width)
            qw = scaled_width(q_block_end - q_block_start + 1, max_query_length, contig_width)
            q_label = "q:{}-{}".format(hit.query_min, hit.query_max)
            q_cx = _clamp_center_x(qx + qw / 2.0, q_label, 9.5, contig_x0, contig_right_limit)
            half_qw = rc.text_width(q_label, 9.5) / 2.0
            q_entries.append(
                {
                    "hit": hit,
                    "x": qx,
                    "w": qw,
                    "label": q_label,
                    "label_cx": q_cx,
                    "span_left": q_cx - half_qw - 2.0,
                    "span_right": q_cx + half_qw + 2.0,
                }
            )

        for hit in sorted(hits, key=lambda h: (h.subject_min, h.subject_max, -h.bitscore)):
            sx = scale(hit.subject_min, protein_length, protein_x0, protein_width)
            sw = scaled_width(hit.subject_max - hit.subject_min + 1, protein_length, protein_width)
            s_label = "s:{}-{}".format(hit.subject_min, hit.subject_max)
            s_cx = _clamp_center_x(sx + sw / 2.0, s_label, 9.5, protein_x0 - 8, protein_right_limit)
            half_sw = rc.text_width(s_label, 9.5) / 2.0
            s_entries.append(
                {
                    "hit": hit,
                    "x": sx,
                    "w": sw,
                    "label": s_label,
                    "label_cx": s_cx,
                    "span_left": s_cx - half_sw - 2.0,
                    "span_right": s_cx + half_sw + 2.0,
                }
            )

        q_lanes = rc.assign_text_lanes(q_entries, pad_px=10.0)
        s_lanes = rc.assign_text_lanes(s_entries, pad_px=10.0)
        left_block_h = len(id_lines) * 13 + len(sub_lines) * 12 + 8
        track_block_h = 34 + max(1, q_lanes, s_lanes) * 15
        row_h = max(62, left_block_h + 14, track_block_h + 12)
        row_layouts.append(
            {
                "qseqid": qseqid,
                "qlen": qlen,
                "anchor": anchor,
                "id_lines": id_lines,
                "sub_lines": sub_lines,
                "q_entries": q_entries,
                "s_entries": s_entries,
                "row_top": cursor_y,
                "row_h": row_h,
            }
        )
        cursor_y += row_h

    total_height = cursor_y + 28
    fig = rc.Figure(fig_width, total_height)

    # Header
    ty = 26
    for line in title_lines:
        fig.text(label_col_x, ty, line, **title_style)
        ty += 17
    ty += 2
    for line in meta_lines:
        fig.text(label_col_x, ty, line, **subtle_style)
        ty += 14

    legend_w = 230
    legend_x = fig_width - 60 - legend_w
    legend_y = 32
    legend_h = 11
    fig.text(
        legend_x + legend_w / 2.0,
        legend_y - 8,
        "Amino-acid identity (pident)",
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

    # Panel a & Panel b Headers + Top Reference / Rulers
    fig.text(contig_x0, panel_header_y, "a", **panel_tag_style)
    fig.text(contig_x0 + 20, panel_header_y, "Contig nucleotide coordinates (nt)", **panel_title_style)

    fig.text(protein_x0, panel_header_y, "b", **panel_tag_style)
    fig.text(protein_x0 + 20, panel_header_y, "Viral protein reference coordinates (aa)", **panel_title_style)

    fig.text(label_col_x, protein_track_y - 1, "Protein reference", **label_style)
    fig.text(label_col_x, protein_track_y + 13, "{} aa reference sequence".format(protein_length), **subtle_style)

    # Panel a shared coordinate ruler (1 .. max_query_length nt)
    fig.polyline([(contig_x0, ruler_y), (contig_x0 + contig_width, ruler_y)], stroke="#555555", stroke_width=1.0)
    last_q_tick_right = -1e9
    for tick_pos in _axis_ticks(max_query_length):
        tx = scale(tick_pos, max_query_length, contig_x0, contig_width)
        fig.polyline([(tx, ruler_y), (tx, ruler_y + 4.0)], stroke="#555555", stroke_width=1.0)
        t_label = "{} nt".format(tick_pos) if tick_pos in (1, max_query_length) else str(tick_pos)
        t_cx = _clamp_center_x(tx, t_label, 9.5, contig_x0, contig_right_limit)
        half_tw = rc.text_width(t_label, 9.5) / 2.0
        if t_cx - half_tw > last_q_tick_right + 10.0:
            fig.text(t_cx, ruler_y + 17.0, t_label, anchor="middle", **axis_style)
            last_q_tick_right = t_cx + half_tw

    # Panel b reference protein bar + shared coordinate ruler (1 .. protein_length aa)
    fig.rect(
        protein_x0,
        protein_track_y - baseline_height / 2.0,
        protein_width,
        baseline_height,
        rx=2,
        fill="#272727",
        stroke=None,
    )
    fig.polyline([(protein_x0, ruler_y), (protein_x0 + protein_width, ruler_y)], stroke="#555555", stroke_width=1.0)
    last_p_tick_right = -1e9
    for tick_pos in _axis_ticks(protein_length):
        tx = scale(tick_pos, protein_length, protein_x0, protein_width)
        fig.polyline([(tx, ruler_y), (tx, ruler_y + 4.0)], stroke="#555555", stroke_width=1.0)
        t_label = "{} aa".format(tick_pos) if tick_pos in (1, protein_length) else str(tick_pos)
        t_cx = _clamp_center_x(tx, t_label, 9.5, protein_x0 - 8, protein_right_limit)
        half_tw = rc.text_width(t_label, 9.5) / 2.0
        if t_cx - half_tw > last_p_tick_right + 10.0:
            fig.text(t_cx, ruler_y + 17.0, t_label, anchor="middle", **axis_style)
            last_p_tick_right = t_cx + half_tw

    # Per-contig alignment rows
    for layout in row_layouts:
        qlen = layout["qlen"]
        anchor = layout["anchor"]
        row_top = layout["row_top"]
        track_y = row_top + 18.0

        # Left column text
        ly = row_top + 10.0
        for line in layout["id_lines"]:
            fig.text(label_col_x, ly, line, **label_style)
            ly += 13.0
        for line in layout["sub_lines"]:
            fig.text(label_col_x, ly, line, **subtle_style)
            ly += 12.0

        # Panel a contig scaffold bar
        contig_bar_w = scaled_width(qlen, max_query_length, contig_width)
        contig_end_x = contig_x0 + contig_bar_w
        fig.rect(
            contig_x0,
            track_y - baseline_height / 2.0,
            contig_bar_w,
            baseline_height,
            rx=2,
            fill="#DCE1E7",
            stroke=None,
        )

        left_coord = qlen if anchor.strand == "-" else 1
        right_coord = 1 if anchor.strand == "-" else qlen
        left_str = "{} nt ({})".format(left_coord, anchor.frame_label)
        right_str = "{} nt".format(right_coord)
        left_w = rc.text_width(left_str, 9.5)
        right_w = rc.text_width(right_str, 9.5)
        rx_end = max(contig_x0 + right_w, min(contig_right_limit, contig_end_x))
        if contig_bar_w >= left_w + right_w + 18.0 and (rx_end - right_w) > (contig_x0 + left_w + 10.0):
            fig.text(contig_x0, track_y - 10.0, left_str, anchor="start", **axis_style)
            fig.text(rx_end, track_y - 10.0, right_str, anchor="end", **axis_style)
        else:
            combined = "{}..{} nt ({})".format(left_coord, right_coord, anchor.frame_label)
            c_cx = _clamp_center_x(contig_x0 + contig_bar_w / 2.0, combined, 9.5, contig_x0, contig_right_limit)
            fig.text(c_cx, track_y - 10.0, combined, anchor="middle", **axis_style)

        # Panel b protein scaffold bar
        fig.rect(
            protein_x0,
            track_y - baseline_height / 2.0,
            protein_width,
            baseline_height,
            rx=2,
            fill="#EAEEF2",
            stroke=None,
        )

        # Panel a HSP blocks & lane-staggered q:start-end labels
        for q_entry in layout["q_entries"]:
            hit = q_entry["hit"]
            color = pident_to_color(hit.pident, identity_low, identity_high)
            fig.rect(
                q_entry["x"],
                track_y - hit_height / 2.0,
                q_entry["w"],
                hit_height,
                rx=2,
                fill=color,
                stroke=None,
                title=hit_title(hit),
            )
            label_y = track_y + baseline_height / 2.0 + 13.0 + q_entry["lane"] * 14.0
            fig.text(q_entry["label_cx"], label_y, q_entry["label"], **coord_style)

        # Panel b HSP blocks & lane-staggered s:start-end labels
        for s_entry in layout["s_entries"]:
            hit = s_entry["hit"]
            color = pident_to_color(hit.pident, identity_low, identity_high)
            fig.rect(
                s_entry["x"],
                track_y - hit_height / 2.0,
                s_entry["w"],
                hit_height,
                rx=2,
                fill=color,
                stroke=None,
                title=hit_title(hit),
            )
            label_y = track_y + baseline_height / 2.0 + 13.0 + s_entry["lane"] * 14.0
            fig.text(s_entry["label_cx"], label_y, s_entry["label"], **coord_style)

    return rc.render_html(fig, protein_label + " blastx summary", protein_label, meta_text)


def write_output_bundle(out_path, html_doc, write_svg=True, write_pdf=True):
    """Delegate to the shared rendering layer (see render_common)."""
    return rc.write_output_bundle(out_path, html_doc, write_svg=write_svg, write_pdf=write_pdf)


def main():
    args = parse_args()
    hits = load_blastx_hits(args.blastx, args.protein, args.min_align_aa)
    if not hits:
        raise SystemExit("No blastx hits found for protein {}".format(args.protein))

    allowed = load_allowed_contigs(args.final_tsv, args.protein)
    if allowed is not None:
        hits = [hit for hit in hits if hit.qseqid in allowed]
        if not hits:
            raise SystemExit("No blastx hits remain after filtering by final candidates for {}".format(args.protein))

    grouped = group_hits_by_contig(hits)
    sorted_names = sorted(
        grouped,
        key=lambda name: (anchor_hit(grouped[name]).subject_min, best_bitscore(grouped[name]), max(hit.qlen for hit in grouped[name])),
        reverse=False,
    )
    if args.top_contigs > 0:
        sorted_names = sorted_names[: args.top_contigs]
    grouped = {name: grouped[name] for name in sorted_names}

    protein_label = load_protein_label(args.protein, args.protein_fasta, args.final_tsv)
    html_doc = build_html(args.protein, protein_label, grouped)
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
