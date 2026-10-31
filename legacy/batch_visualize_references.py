#!/usr/bin/env python3

import argparse
import csv
import importlib.util
import re
from pathlib import Path


def load_visualizer_module(module_path):
    spec = importlib.util.spec_from_file_location("visualize_blastn_reference", str(module_path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch-generate blastn reference graphic summaries (HTML + SVG + PDF) from final_retained_candidates.tsv."
    )
    parser.add_argument(
        "--final-tsv",
        default="final_retained_candidates.tsv",
        help="Input final_retained_candidates.tsv",
    )
    parser.add_argument(
        "--blastn",
        default="blastn_res.txt",
        help="Input blastn result file",
    )
    parser.add_argument(
        "--outdir",
        default="reference_graphic_summaries",
        help="Output directory for generated HTML/PDF files",
    )
    parser.add_argument(
        "--module",
        default="visualize_blastn_reference.py",
        help="Path to visualize_blastn_reference.py",
    )
    parser.add_argument(
        "--top-contigs",
        type=int,
        default=0,
        help="Optional limit on contigs displayed in each figure (0 means all).",
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
        help="Skip PDF export.",
    )
    parser.add_argument(
        "--no-svg",
        action="store_true",
        help="Skip separate SVG export.",
    )
    return parser.parse_args()


def extract_references(final_tsv):
    references = []
    seen = set()
    with Path(final_tsv).open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            best_nt_hit = (row.get("best_nt_hit") or "").strip()
            if not best_nt_hit or best_nt_hit.startswith("no strong nucleotide hit"):
                continue
            ref = best_nt_hit.split(" | ", 1)[0].strip()
            if ref and ref not in seen:
                seen.add(ref)
                references.append(ref)
    return references


def limit_grouped_hits(module, grouped_hits, top_contigs):
    if top_contigs <= 0:
        return grouped_hits
    sorted_names = sorted(
        grouped_hits,
        key=lambda name: (
            module.best_bitscore(grouped_hits[name]),
            max(hit.qlen for hit in grouped_hits[name]),
        ),
        reverse=True,
    )
    sorted_names = sorted_names[:top_contigs]
    return {name: grouped_hits[name] for name in sorted_names}


def reference_filename_prefix(reference_label, reference):
    if " | " in reference_label:
        name = reference_label.split(" | ", 1)[1]
    else:
        name = reference_label
    name = re.split(r",\s*(?:complete|coding-complete|partial)\b", name, maxsplit=1)[0]
    slug = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
    if not slug:
        slug = reference
    return "{}__{}".format(reference, slug)


def write_manifest(path, rows):
    fieldnames = [
        "reference",
        "reference_label",
        "html_path",
        "svg_path",
        "pdf_path",
        "contig_count",
        "hsp_count",
        "status",
    ]
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    module = load_visualizer_module(Path(args.module).resolve())
    references = extract_references(args.final_tsv)
    manifest_rows = []

    for reference in references:
        hits = module.load_blastn_hits(args.blastn, reference, args.min_align_len)
        allowed = module.load_allowed_contigs(args.final_tsv, reference)
        hits = [hit for hit in hits if hit.qseqid in allowed]

        if not hits:
            manifest_rows.append(
                {
                    "reference": reference,
                    "reference_label": reference,
                    "html_path": "",
                    "svg_path": "",
                    "pdf_path": "",
                    "contig_count": "0",
                    "hsp_count": "0",
                    "status": "skipped_no_hits",
                }
            )
            continue

        grouped_hits = module.group_hits_by_contig(hits)
        grouped_hits = limit_grouped_hits(module, grouped_hits, args.top_contigs)
        reference_label = module.load_reference_label(args.final_tsv, reference)

        prefix = reference_filename_prefix(reference_label, reference)
        html_path = outdir / (prefix + ".graphic_summary.html")
        html_text = module.build_html(reference, reference_label, grouped_hits)
        html_path, svg_path, pdf_path = module.write_output_bundle(
            html_path,
            html_text,
            write_svg=not args.no_svg,
            write_pdf=not args.no_pdf,
        )

        manifest_rows.append(
            {
                "reference": reference,
                "reference_label": reference_label,
                "html_path": str(html_path.resolve()),
                "svg_path": str(svg_path.resolve()) if svg_path is not None else "",
                "pdf_path": str(pdf_path.resolve()) if pdf_path is not None else "",
                "contig_count": str(len(grouped_hits)),
                "hsp_count": str(sum(len(hits) for hits in grouped_hits.values())),
                "status": "ok",
            }
        )

        outputs = [html_path.name]
        if svg_path is not None:
            outputs.append(svg_path.name)
        if pdf_path is not None:
            outputs.append(pdf_path.name)
        print("Generated:", reference, "->", ", ".join(outputs))

    manifest_path = outdir / "reference_graphic_summary_manifest.tsv"
    write_manifest(manifest_path, manifest_rows)
    print("Manifest:", manifest_path.resolve())


if __name__ == "__main__":
    main()
