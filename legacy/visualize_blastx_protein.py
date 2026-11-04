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


def build_html(protein, protein_label, grouped_hits):
    """Build the two-panel figure (contig nt / protein aa) via the shared
    rendering layer. Panel geometry and frame annotations are specific to
    the protein view; drawing primitives come from rc.Figure so the same
    figure renders to inline SVG (HTML) and native vector PDF.
    """
    protein_length = next(iter(next(iter(grouped_hits.values())))).slen if grouped_hits else 0
    max_query_length = max((max(hit.qlen for hit in hits) for hits in grouped_hits.values()), default=0)

    sorted_items = sorted(
        grouped_hits.items(),
        key=lambda item: (anchor_hit(item[1]).subject_min, -best_bitscore(item[1]), -max(hit.qlen for hit in item[1])),
    )

    identity_low, identity_high = observed_identity_range(grouped_hits)

    label_style = {"size": 12, "fill": "#222222"}
    subtle_style = {"size": 11, "fill": "#666666"}
    axis_style = {"size": 11, "fill": "#555555"}
    legend_style = {"size": 11, "fill": "#444444"}
    coord_style = {"size": 9, "fill": "#333333", "anchor": "middle"}

    contig_x0 = 280
    contig_width = 620
    protein_x0 = 960
    protein_width = 360
    row_height = 64
    protein_track_y = 92
    baseline_height = 14
    hit_height = 10
    header_height = 150
    total_height = header_height + len(sorted_items) * row_height + 40

    fig = rc.Figure(1400, total_height)

    meta_text = (
        "Protein length: {} aa | Longest contig: {} nt | Displayed contigs: {} | Displayed HSPs: {} | "
        "pident range: {:.3f}-{:.3f}"
    ).format(
        protein_length,
        max_query_length,
        len(sorted_items),
        sum(len(hits) for _, hits in sorted_items),
        identity_low,
        identity_high,
    )

    fig.text(24, 28, "Left panel: contig nt coordinates. Right panel: protein aa coordinates.", **subtle_style)

    fig.gradient_rect(
        1080,
        30,
        220,
        12,
        lambda frac: pident_to_color(identity_low + (identity_high - identity_low) * frac, identity_low, identity_high),
        stroke="#bbbbbb",
    )
    fig.text(1080, 58, "{:.2f}%".format(identity_low), **legend_style)
    fig.text(1264, 58, "{:.2f}%".format(identity_high), **legend_style)
    fig.text(1172, 58, "pident", **legend_style)

    fig.text(24, protein_track_y + 4, "Protein reference", **label_style)
    fig.rect(
        protein_x0,
        protein_track_y - baseline_height / 2,
        protein_width,
        baseline_height,
        rx=2,
        fill="#111111",
        stroke="#111111",
        stroke_width=1.0,
    )
    fig.text(protein_x0 - 10, protein_track_y + 24, "1 aa", **axis_style)
    fig.text(protein_x0 + protein_width - 42, protein_track_y + 24, "{} aa".format(protein_length), **axis_style)

    current_y = protein_track_y + 56
    for qseqid, hits in sorted_items:
        qlen = max(hit.qlen for hit in hits)
        best = best_bitscore(hits)
        anchor = anchor_hit(hits)
        fig.text(24, current_y + 4, qseqid, **label_style)
        fig.text(
            24,
            current_y + 18,
            "len={} nt | HSPs={} | best bitscore={:.1f} | strand={} | frame={}".format(
                qlen,
                len(hits),
                best,
                anchor.strand,
                anchor.frame_label,
            ),
            **subtle_style
        )

        contig_y = current_y
        protein_y = current_y

        fig.rect(
            contig_x0,
            contig_y - baseline_height / 2,
            scaled_width(qlen, max_query_length, contig_width),
            baseline_height,
            rx=2,
            fill="#a0a0a0",
            stroke="#8d8d8d",
            stroke_width=1.0,
        )
        contig_end_x = contig_x0 + scaled_width(qlen, max_query_length, contig_width)
        left_label = "{} nt".format(qlen) if anchor.strand == "-" else "1 nt"
        right_label = "1 nt" if anchor.strand == "-" else "{} nt".format(qlen)
        fig.text(contig_x0 - 22, contig_y + 24, left_label, **axis_style)
        fig.text(contig_end_x - 36, contig_y + 24, right_label, **axis_style)

        fig.rect(
            protein_x0,
            protein_y - baseline_height / 2,
            protein_width,
            baseline_height,
            rx=2,
            fill="#d9d9d9",
            stroke="#b0b0b0",
            stroke_width=1.0,
        )
        fig.text(protein_x0 - 10, protein_y + 24, "1 aa", **axis_style)
        fig.text(protein_x0 + protein_width - 42, protein_y + 24, "{} aa".format(protein_length), **axis_style)

        for hit in hits:
            color = pident_to_color(hit.pident, identity_low, identity_high)
            title = hit_title(hit)

            if anchor.strand == "+":
                q_block_start = hit.query_min
                q_block_end = hit.query_max
            else:
                q_block_start = qlen - hit.query_max + 1
                q_block_end = qlen - hit.query_min + 1

            qx = scale(q_block_start, max_query_length, contig_x0, contig_width)
            qw = scaled_width(q_block_end - q_block_start + 1, max_query_length, contig_width)
            sx = scale(hit.subject_min, protein_length, protein_x0, protein_width)
            sw = scaled_width(hit.subject_max - hit.subject_min + 1, protein_length, protein_width)

            fig.rect(
                qx,
                contig_y - hit_height / 2,
                qw,
                hit_height,
                rx=2,
                fill=color,
                title=title,
            )
            fig.rect(
                sx,
                protein_y - hit_height / 2,
                sw,
                hit_height,
                rx=2,
                fill=color,
                title=title,
            )

            fig.text(
                qx + qw / 2.0,
                contig_y + baseline_height / 2 + 12,
                "q:{}-{}".format(hit.query_min, hit.query_max),
                **coord_style
            )
            fig.text(
                sx + sw / 2.0,
                protein_y - baseline_height / 2 - 6,
                "s:{}-{}".format(hit.subject_min, hit.subject_max),
                **coord_style
            )

        current_y += row_height

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
