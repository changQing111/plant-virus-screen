#!/usr/bin/env python3

import argparse
import csv
import html
import math
import re
import shutil
import subprocess
import colorsys
from pathlib import Path


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
        description="Visualize blastn hits for one reference accession as a self-contained HTML/SVG summary."
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
        help="Only write HTML (and optional PDF) and skip separate SVG export.",
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
    # Low identity stays blue/cyan, high identity shifts to red.
    hue = (220.0 - 220.0 * ratio) / 360.0
    saturation = 0.78
    lightness = 0.50
    r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
    return "#{:02X}{:02X}{:02X}".format(int(round(r * 255)), int(round(g * 255)), int(round(b * 255)))


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
    reference_length = next(iter(next(iter(grouped_hits.values())))).slen if grouped_hits else 0
    left_margin = 220
    right_margin = 80
    track_width = 1100
    row_height = 54
    baseline_height = 14
    hit_height = baseline_height
    gap_after_reference = 42
    lane_step = 14

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

    parts = []
    parts.append("<!DOCTYPE html>")
    parts.append("<html lang='en'>")
    parts.append("<head>")
    parts.append("<meta charset='utf-8'>")
    parts.append("<title>{}</title>".format(html.escape(reference_label + " blastn summary")))
    parts.append(
        "<style>"
        "body{font-family:Arial,Helvetica,sans-serif;margin:24px;color:#222;background:#fff;}"
        "h1{font-size:22px;margin:0 0 8px 0;}"
        ".meta{margin:0 0 18px 0;font-size:14px;color:#444;}"
        "svg{border:1px solid #d0d0d0;background:#fafafa;}"
        ".label{font-size:12px;fill:#222;}"
        ".subtle{font-size:11px;fill:#666;}"
        ".axis{font-size:11px;fill:#555;}"
        ".baseline-ref{fill:#111;stroke:#111;stroke-width:1;}"
        ".baseline-query{fill:#a0a0a0;stroke:#8d8d8d;stroke-width:1;}"
        ".ref-marker{stroke:#111;stroke-width:1.5;fill:none;}"
        ".legend-text{font-size:11px;fill:#444;}"
        ".coord-text{font-size:9px;fill:#333;}"
        "</style>"
    )
    parts.append("</head>")
    parts.append("<body>")
    parts.append("<h1>{}</h1>".format(html.escape(reference_label)))
    parts.append(
        "<div class='meta'>Reference length: {} bp | Longest contig: {} bp | Displayed contigs: {} | Displayed HSPs: {} | "
        "Display axis: {}..{} bp | "
        "pident range: {:.3f}-{:.3f}</div>".format(
            reference_length,
            max_query_length,
            len(sorted_items),
            sum(len(hits) for _, hits in sorted_items),
            axis_min,
            axis_max,
            identity_low,
            identity_high,
        )
    )

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

    parts.append("<svg width='{}' height='{}' viewBox='0 0 {} {}'>".format(left_margin + track_width + right_margin, total_height, left_margin + track_width + right_margin, total_height))

    parts.append("<text x='20' y='28' class='subtle'>Reference track in black, contig tracks in gray. Colored blocks mark HSPs.</text>")

    legend_x = left_margin + track_width - 260
    legend_y = 30
    legend_w = 220
    legend_h = 12
    gradient_id = "pid-gradient"
    gradient_steps = 64
    parts.append("<defs>")
    parts.append("<linearGradient id='{}' x1='0%' y1='0%' x2='100%' y2='0%'>".format(gradient_id))
    for idx in range(gradient_steps):
        frac = float(idx) / max(1, gradient_steps - 1)
        pident = identity_low + (identity_high - identity_low) * frac
        color = pident_to_color(pident, identity_low, identity_high)
        parts.append("<stop offset='{:.2f}%' stop-color='{}'/>".format(frac * 100.0, color))
    parts.append("</linearGradient>")
    parts.append("</defs>")
    parts.append("<rect x='{:.2f}' y='{}' width='{:.2f}' height='{}' fill='url(#{})' stroke='#bbb'/>".format(
        legend_x,
        legend_y,
        legend_w,
        legend_h,
        gradient_id,
    ))
    parts.append("<text x='{:.2f}' y='{}' class='legend-text'>{:.2f}%</text>".format(legend_x, legend_y + 28, identity_low))
    parts.append("<text x='{:.2f}' y='{}' class='legend-text'>{:.2f}%</text>".format(legend_x + legend_w - 36, legend_y + 28, identity_high))
    parts.append("<text x='{:.2f}' y='{}' class='legend-text'>pident</text>".format(legend_x + legend_w / 2 - 18, legend_y + 28))

    ref_marker_y = ref_y - baseline_height / 2 - 7
    parts.append("<text x='20' y='{}' class='label'>{}</text>".format(ref_y + 4, html.escape(reference)))
    ref_start_x = scaled_x(1, axis_min, axis_span, ref_x0, track_width)
    ref_bar_width = scaled_width(reference_length, axis_span, track_width)
    parts.append(
        "<rect x='{:.2f}' y='{}' width='{:.2f}' height='{}' rx='2' class='baseline-ref'/>".format(
            ref_start_x,
            ref_y - baseline_height / 2,
            ref_bar_width,
            baseline_height,
        )
    )
    parts.append("<text x='{:.2f}' y='{}' class='axis'>1</text>".format(ref_start_x - 8, ref_y + 24))
    ref_end_x = ref_start_x + ref_bar_width
    parts.append("<text x='{:.2f}' y='{}' class='axis'>{}</text>".format(ref_end_x - 30, ref_y + 24, reference_length))

    # Keep the reference bar pure black, but mark each HSP's true interval with
    # a black bracket above the bar so the subject coordinates remain visible.
    for entry in reference_entries:
        hit = entry["hit"]
        ref_hit_x1 = entry["x1"]
        ref_hit_x2 = entry["x2"]
        title = html.escape(title_text(hit))
        ref_mid_x = ref_hit_x1 + entry["width"] / 2.0
        lane = entry["lane"]
        lane_marker_y = ref_marker_y - lane * lane_step
        parts.append(
            "<path d='M {x1:.2f} {y1:.2f} V {y2:.2f} H {x2:.2f} V {y1:.2f}' class='ref-marker'>"
            "<title>{title}</title></path>".format(
                x1=ref_hit_x1,
                y1=lane_marker_y,
                y2=ref_y - baseline_height / 2,
                x2=ref_hit_x2,
                title=title,
            )
        )
        parts.append(
            "<text x='{:.2f}' y='{:.2f}' text-anchor='middle' class='coord-text'>s:{}-{}</text>".format(
                ref_mid_x,
                lane_marker_y - 3,
                hit.subject_min,
                hit.subject_max,
            )
        )

    current_y = ref_y + gap_after_reference
    for qseqid, hits in sorted_items:
        qlen = contig_positions[qseqid]["qlen"]
        contig_start_abs = contig_positions[qseqid]["contig_start_abs"]
        contig_end_abs = contig_positions[qseqid]["contig_end_abs"]
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
        parts.append("<text x='20' y='{}' class='label'>{}</text>".format(current_y + 4, html.escape(label_main)))
        parts.append("<text x='20' y='{}' class='subtle'>{}</text>".format(current_y + 18, html.escape(label_sub)))

        query_y = current_y
        query_start_x = scaled_x(contig_start_abs, axis_min, axis_span, ref_x0, track_width)
        query_bar_width = scaled_width(qlen, axis_span, track_width)
        parts.append(
            "<rect x='{:.2f}' y='{}' width='{:.2f}' height='{}' rx='2' class='baseline-query'/>".format(
                query_start_x,
                query_y - baseline_height / 2,
                query_bar_width,
                baseline_height,
            )
        )
        left_label = qlen if contig_orientation == "-" else 1
        right_label = 1 if contig_orientation == "-" else qlen
        parts.append("<text x='{:.2f}' y='{}' class='axis'>{}</text>".format(query_start_x - 16, query_y + 24, left_label))
        query_end_x = query_start_x + query_bar_width
        parts.append("<text x='{:.2f}' y='{}' class='axis'>{}</text>".format(query_end_x - 30, query_y + 24, right_label))

        for hit in hits:
            color = pident_to_color(hit.pident, identity_low, identity_high)
            qry_hit_start_abs, qry_hit_end_abs = query_interval_absolute(
                hit, contig_start_abs, qlen, contig_orientation
            )
            qry_hit_x1 = scaled_x(qry_hit_start_abs, axis_min, axis_span, ref_x0, track_width)
            title = html.escape(title_text(hit))
            qry_hit_width = scaled_width(qry_hit_end_abs - qry_hit_start_abs + 1, axis_span, track_width)
            qry_mid_x = qry_hit_x1 + qry_hit_width / 2.0

            parts.append(
                "<rect x='{:.2f}' y='{}' width='{:.2f}' height='{}' rx='2' fill='{}'>"
                "<title>{}</title></rect>".format(
                    qry_hit_x1,
                    query_y - hit_height / 2,
                    qry_hit_width,
                    hit_height,
                    color,
                    title,
                )
            )
            parts.append(
                "<text x='{:.2f}' y='{:.2f}' text-anchor='middle' class='coord-text'>q:{}-{}</text>".format(
                    qry_mid_x,
                    query_y + baseline_height / 2 + 12,
                    hit.query_min,
                    hit.query_max,
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
    outdir = str(html_path.parent)
    cmd = [
        tool,
        "--headless",
        "--convert-to",
        "pdf",
        "--outdir",
        outdir,
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
        raise RuntimeError(
            "PDF export command finished without creating {}.\nCommand output:\n{}".format(
                pdf_path, result.stdout
            )
        )
    return pdf_path


def export_pdf_from_svg(svg_path):
    tool = shutil.which("soffice") or shutil.which("libreoffice")
    if not tool:
        raise RuntimeError("Neither soffice nor libreoffice was found; cannot export PDF automatically.")

    svg_path = Path(svg_path).resolve()
    pdf_path = svg_path.with_suffix(".pdf")
    if pdf_path.exists():
        pdf_path.unlink()
    outdir = str(svg_path.parent)
    cmd = [
        tool,
        "--headless",
        "--convert-to",
        "pdf",
        "--outdir",
        outdir,
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
        raise RuntimeError("PDF export failed:\n{}".format(result.stdout))
    if not pdf_path.exists():
        raise RuntimeError(
            "PDF export command finished without creating {}.\nCommand output:\n{}".format(
                pdf_path, result.stdout
            )
        )
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
    html_text = build_html(args.reference, reference_label, grouped)
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
