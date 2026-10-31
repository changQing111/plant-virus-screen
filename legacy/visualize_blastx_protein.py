#!/usr/bin/env python3

import argparse
import csv
import html
import math
import os
import re
import shutil
import subprocess
import colorsys
from pathlib import Path


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
        description="Visualize blastx hits for one viral protein accession as a self-contained HTML summary with optional PDF export."
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
        help="Only write HTML (and optional PDF) and skip separate SVG export.",
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
    values = [hit.pident for hits in grouped_hits.values() for hit in hits]
    if not values:
        return 0.0, 100.0
    low = min(values)
    high = max(values)
    if math.isclose(low, high):
        low = max(0.0, low - 5.0)
        high = min(100.0, high + 5.0)
    return low, high


def pident_to_color(pident, low, high):
    if high <= low:
        ratio = 0.5
    else:
        ratio = (pident - low) / (high - low)
    ratio = max(0.0, min(1.0, ratio))
    hue = (220.0 - 220.0 * ratio) / 360.0
    saturation = 0.78
    lightness = 0.50
    r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
    return "#{:02X}{:02X}{:02X}".format(int(round(r * 255)), int(round(g * 255)), int(round(b * 255)))


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
    protein_length = next(iter(next(iter(grouped_hits.values())))).slen if grouped_hits else 0
    max_query_length = max((max(hit.qlen for hit in hits) for hits in grouped_hits.values()), default=0)

    sorted_items = sorted(
        grouped_hits.items(),
        key=lambda item: (anchor_hit(item[1]).subject_min, -best_bitscore(item[1]), -max(hit.qlen for hit in item[1])),
    )

    identity_low, identity_high = observed_identity_range(grouped_hits)

    left_margin = 24
    label_width = 250
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

    parts = []
    parts.append("<!DOCTYPE html>")
    parts.append("<html lang='en'>")
    parts.append("<head>")
    parts.append("<meta charset='utf-8'>")
    parts.append("<title>{}</title>".format(html.escape(protein_label + " blastx summary")))
    parts.append(
        "<style>"
        "body{font-family:Arial,Helvetica,sans-serif;margin:24px;color:#222;background:#fff;}"
        "h1{font-size:22px;margin:0 0 8px 0;}"
        ".meta{margin:0 0 18px 0;font-size:14px;color:#444;}"
        "svg{border:1px solid #d0d0d0;background:#fafafa;}"
        ".label{font-size:12px;fill:#222;}"
        ".subtle{font-size:11px;fill:#666;}"
        ".axis{font-size:11px;fill:#555;}"
        ".protein-ref{fill:#111;stroke:#111;stroke-width:1;}"
        ".protein-track{fill:#d9d9d9;stroke:#b0b0b0;stroke-width:1;}"
        ".contig-track{fill:#a0a0a0;stroke:#8d8d8d;stroke-width:1;}"
        ".legend-text{font-size:11px;fill:#444;}"
        ".coord-text{font-size:9px;fill:#333;}"
        "</style>"
    )
    parts.append("</head>")
    parts.append("<body>")
    parts.append("<h1>{}</h1>".format(html.escape(protein_label)))
    parts.append(
        "<div class='meta'>Protein length: {} aa | Longest contig: {} nt | Displayed contigs: {} | Displayed HSPs: {} | "
        "pident range: {:.3f}-{:.3f}</div>".format(
            protein_length,
            max_query_length,
            len(sorted_items),
            sum(len(hits) for _, hits in sorted_items),
            identity_low,
            identity_high,
        )
    )
    parts.append("<svg width='1400' height='{}' viewBox='0 0 1400 {}'>".format(total_height, total_height))
    parts.append("<text x='24' y='28' class='subtle'>Left panel: contig nt coordinates. Right panel: protein aa coordinates.</text>")

    gradient_id = "pid-gradient"
    parts.append("<defs>")
    parts.append("<linearGradient id='{}' x1='0%' y1='0%' x2='100%' y2='0%'>".format(gradient_id))
    for pct in (0, 25, 50, 75, 100):
        pident = identity_low + (identity_high - identity_low) * (pct / 100.0)
        parts.append("<stop offset='{}%' stop-color='{}'/>".format(pct, pident_to_color(pident, identity_low, identity_high)))
    parts.append("</linearGradient>")
    parts.append("</defs>")

    parts.append("<rect x='1080' y='30' width='220' height='12' fill='url(#{})' stroke='#bbb'/>".format(gradient_id))
    parts.append("<text x='1080' y='58' class='legend-text'>{:.2f}%</text>".format(identity_low))
    parts.append("<text x='1264' y='58' class='legend-text'>{:.2f}%</text>".format(identity_high))
    parts.append("<text x='1172' y='58' class='legend-text'>pident</text>")

    parts.append("<text x='24' y='{}' class='label'>Protein reference</text>".format(protein_track_y + 4))
    parts.append("<rect x='{:.2f}' y='{}' width='{:.2f}' height='{}' rx='2' class='protein-ref'/>".format(
        protein_x0,
        protein_track_y - baseline_height / 2,
        protein_width,
        baseline_height,
    ))
    parts.append("<text x='{:.2f}' y='{}' class='axis'>1 aa</text>".format(protein_x0 - 10, protein_track_y + 24))
    parts.append("<text x='{:.2f}' y='{}' class='axis'>{} aa</text>".format(protein_x0 + protein_width - 42, protein_track_y + 24, protein_length))

    current_y = protein_track_y + 56
    for qseqid, hits in sorted_items:
        qlen = max(hit.qlen for hit in hits)
        best = best_bitscore(hits)
        anchor = anchor_hit(hits)
        parts.append("<text x='24' y='{}' class='label'>{}</text>".format(current_y + 4, html.escape(qseqid)))
        parts.append(
            "<text x='24' y='{}' class='subtle'>len={} nt | HSPs={} | best bitscore={:.1f} | strand={} | frame={}</text>".format(
                current_y + 18,
                qlen,
                len(hits),
                best,
                anchor.strand,
                anchor.frame_label,
            )
        )

        contig_y = current_y
        protein_y = current_y

        parts.append("<rect x='{:.2f}' y='{}' width='{:.2f}' height='{}' rx='2' class='contig-track'/>".format(
            contig_x0,
            contig_y - baseline_height / 2,
            scaled_width(qlen, max_query_length, contig_width),
            baseline_height,
        ))
        contig_end_x = contig_x0 + scaled_width(qlen, max_query_length, contig_width)
        left_label = "{} nt".format(qlen) if anchor.strand == "-" else "1 nt"
        right_label = "1 nt" if anchor.strand == "-" else "{} nt".format(qlen)
        parts.append("<text x='{:.2f}' y='{}' class='axis'>{}</text>".format(contig_x0 - 22, contig_y + 24, left_label))
        parts.append("<text x='{:.2f}' y='{}' class='axis'>{}</text>".format(contig_end_x - 36, contig_y + 24, right_label))

        parts.append("<rect x='{:.2f}' y='{}' width='{:.2f}' height='{}' rx='2' class='protein-track'/>".format(
            protein_x0,
            protein_y - baseline_height / 2,
            protein_width,
            baseline_height,
        ))
        parts.append("<text x='{:.2f}' y='{}' class='axis'>1 aa</text>".format(protein_x0 - 10, protein_y + 24))
        parts.append("<text x='{:.2f}' y='{}' class='axis'>{} aa</text>".format(protein_x0 + protein_width - 42, protein_y + 24, protein_length))

        for hit in hits:
            color = pident_to_color(hit.pident, identity_low, identity_high)
            title = html.escape(hit_title(hit))

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

            parts.append(
                "<rect x='{:.2f}' y='{}' width='{:.2f}' height='{}' rx='2' fill='{}'><title>{}</title></rect>".format(
                    qx,
                    contig_y - hit_height / 2,
                    qw,
                    hit_height,
                    color,
                    title,
                )
            )
            parts.append(
                "<rect x='{:.2f}' y='{}' width='{:.2f}' height='{}' rx='2' fill='{}'><title>{}</title></rect>".format(
                    sx,
                    protein_y - hit_height / 2,
                    sw,
                    hit_height,
                    color,
                    title,
                )
            )

            parts.append(
                "<text x='{:.2f}' y='{:.2f}' text-anchor='middle' class='coord-text'>q:{}-{}</text>".format(
                    qx + qw / 2.0,
                    contig_y + baseline_height / 2 + 12,
                    hit.query_min,
                    hit.query_max,
                )
            )
            parts.append(
                "<text x='{:.2f}' y='{:.2f}' text-anchor='middle' class='coord-text'>s:{}-{}</text>".format(
                    sx + sw / 2.0,
                    protein_y - baseline_height / 2 - 6,
                    hit.subject_min,
                    hit.subject_max,
                )
            )

        current_y += row_height

    parts.append("</svg>")
    parts.append("</body>")
    parts.append("</html>")
    return "\n".join(parts)


def extract_svg_text(html_text):
    svg_match = re.search(r"(<svg\b.*?</svg>)", html_text, flags=re.S)
    if not svg_match:
        raise RuntimeError("Failed to extract <svg> block from generated HTML.")
    svg_text = svg_match.group(1)

    style_match = re.search(r"<style>(.*?)</style>", html_text, flags=re.S)
    if style_match:
        css = style_match.group(1)
        css = re.sub(r"body\{.*?\}", "", css)
        css = re.sub(r"h1\{.*?\}", "", css)
        css = re.sub(r"\.meta\{.*?\}", "", css)
        css = re.sub(r"svg\{.*?\}", "", css)
        css = re.sub(r"\s+", " ", css).strip()
        style_block = "<style type='text/css'>{}</style>".format(css)
        svg_text = re.sub(r"(<svg\b[^>]*>)", r"\1\n" + style_block, svg_text, count=1)

    if "xmlns=" not in svg_text.split(">", 1)[0]:
        svg_text = re.sub(r"<svg\b", "<svg xmlns='http://www.w3.org/2000/svg'", svg_text, count=1)

    return "<?xml version='1.0' encoding='UTF-8'?>\n" + svg_text


def export_pdf_from_html(html_path):
    tool = shutil.which("soffice") or shutil.which("libreoffice")
    if not tool:
        raise RuntimeError("Neither soffice nor libreoffice was found; cannot export PDF automatically.")

    html_path = Path(html_path).resolve()
    pdf_path = html_path.with_suffix(".pdf")
    if pdf_path.exists():
        pdf_path.unlink()

    cmd = [
        tool,
        "--headless",
        "--convert-to",
        "pdf",
        "--outdir",
        str(html_path.parent),
        str(html_path),
    ]
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError("PDF export failed:\n{}".format(result.stdout))
    if not pdf_path.exists():
        raise RuntimeError("PDF export finished without creating {}".format(pdf_path))
    return pdf_path


def export_pdf_from_svg(svg_path):
    tool = shutil.which("soffice") or shutil.which("libreoffice")
    if not tool:
        raise RuntimeError("Neither soffice nor libreoffice was found; cannot export PDF automatically.")

    svg_path = Path(svg_path).resolve()
    pdf_path = svg_path.with_suffix(".pdf")
    if pdf_path.exists():
        pdf_path.unlink()

    cmd = [
        tool,
        "--headless",
        "--convert-to",
        "pdf",
        "--outdir",
        str(svg_path.parent),
        str(svg_path),
    ]
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError("PDF export from SVG failed:\n{}".format(result.stdout))
    if not pdf_path.exists():
        raise RuntimeError("PDF export finished without creating {}".format(pdf_path))
    return pdf_path


def write_output_bundle(out_path, html_text, write_svg=True, write_pdf=True):
    out_path = Path(out_path)
    out_path.write_text(html_text, encoding="utf-8")

    svg_path = None
    if write_svg:
        svg_path = out_path.with_suffix(".svg")
        svg_path.write_text(extract_svg_text(html_text), encoding="utf-8")

    pdf_path = None
    if write_pdf:
        if svg_path is not None:
            pdf_path = export_pdf_from_svg(svg_path)
        else:
            pdf_path = export_pdf_from_html(out_path)

    return out_path, svg_path, pdf_path


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
    html_text = build_html(args.protein, protein_label, grouped)
    out_path, svg_path, pdf_path = write_output_bundle(
        args.out,
        html_text,
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
