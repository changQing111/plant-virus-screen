#!/usr/bin/env python3
"""Batch-generate candidate contig coverage figures (HTML + native vector PDF).

Reads the report plot-candidates TSV, scans the samtools depth file once
for all requested contigs, and writes one MAPQ>=20 coverage curve figure
per contig into ``coverage_plots/`` plus a manifest TSV. ``--max-contigs``
caps the number of figures (0 means all) so ``--plots all`` cannot explode
into thousands of files; capped contigs are still listed in the manifest
with status ``skipped_cap`` for auditability.
"""

import argparse
import csv
import importlib.util
import re
from pathlib import Path


def load_visualizer_module(module_path):
    spec = importlib.util.spec_from_file_location("visualize_coverage", str(module_path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch-generate candidate coverage summaries (HTML + native vector PDF) from the report plot-candidates TSV."
    )
    parser.add_argument(
        "--final-tsv",
        default="plot_candidates.tsv",
        help="Input plot candidates TSV (query_id plus coverage/annotation columns)",
    )
    parser.add_argument(
        "--depth",
        default="depth_mapq20.tsv",
        help="samtools depth -a output (contig, position, depth), MAPQ-filtered",
    )
    parser.add_argument(
        "--outdir",
        default="coverage_plots",
        help="Output directory for generated HTML/PDF files",
    )
    parser.add_argument(
        "--module",
        default="visualize_coverage.py",
        help="Path to visualize_coverage.py",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=1000,
        help="Maximum number of curve bins per contig (default 1000).",
    )
    parser.add_argument(
        "--max-contigs",
        type=int,
        default=200,
        help="Cap on the number of contigs plotted, in TSV order (0 means all). Capped contigs are recorded in the manifest.",
    )
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="Skip PDF export.",
    )
    parser.add_argument(
        "--no-svg",
        action="store_true",
        help="Deprecated no-op kept for compatibility: standalone SVG output was replaced by native vector PDF.",
    )
    return parser.parse_args()


def read_candidate_rows(final_tsv):
    rows = []
    with Path(final_tsv).open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            contig = (row.get("query_id") or row.get("contig_id") or "").strip()
            if contig:
                rows.append({"contig_id": contig, "row": row})
    return rows


def sanitize_filename(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def write_manifest(path, rows):
    fieldnames = [
        "contig_id",
        "length",
        "mean_depth",
        "breadth_1x",
        "breadth_5x",
        "evidence_grade",
        "html_path",
        "svg_path",
        "pdf_path",
        "status",
    ]
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def annotation_from_row(row):
    return {
        "evidence_grade": (row.get("evidence_grade") or "").strip(),
        "putative_virus": (row.get("putative_virus") or row.get("candidate_name") or "").strip(),
        "blastn_anno": (row.get("blastn_anno") or "").strip(),
        "blastx_anno": (row.get("blastx_anno") or "").strip(),
        "mean_depth_tsv": (row.get("mapq20_mean_depth") or "").strip(),
        "breadth_1x_tsv": (row.get("mapq20_breadth_1x") or "").strip(),
    }


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    module = load_visualizer_module(Path(args.module).resolve())
    candidates = read_candidate_rows(args.final_tsv)

    kept = candidates
    capped = []
    if args.max_contigs > 0 and len(candidates) > args.max_contigs:
        kept = candidates[: args.max_contigs]
        capped = candidates[args.max_contigs :]

    kept_ids = [item["contig_id"] for item in kept]
    print("Loading depth for {} contigs (single scan of {})".format(len(kept_ids), args.depth))
    depths = module.load_depth_positions(args.depth, kept_ids)

    manifest_rows = []
    generated = 0
    for item in kept:
        contig = item["contig_id"]
        values = depths.get(contig) or []
        base = {
            "contig_id": contig,
            "length": "",
            "mean_depth": "",
            "breadth_1x": "",
            "breadth_5x": "",
            "evidence_grade": (item["row"].get("evidence_grade") or "").strip(),
            "html_path": "",
            "svg_path": "",
            "pdf_path": "",
        }
        if not values:
            manifest_rows.append(dict(base, status="skipped_no_depth"))
            continue

        stats = module.compute_stats(values)
        base.update(
            {
                "length": str(stats["length"]),
                "mean_depth": "{:.4f}".format(stats["mean"]),
                "breadth_1x": "{:.6f}".format(stats["breadth_1x"]),
                "breadth_5x": "{:.6f}".format(stats["breadth_5x"]),
            }
        )

        prefix = sanitize_filename(contig)
        html_path = outdir / (prefix + ".coverage_summary.html")
        annotation = annotation_from_row(item["row"])
        html_doc = module.build_html(contig, values, annotation, bins=args.bins)
        html_path, svg_path, pdf_path = module.write_output_bundle(
            html_path,
            html_doc,
            write_svg=not args.no_svg,
            write_pdf=not args.no_pdf,
        )
        manifest_rows.append(
            dict(
                base,
                html_path=str(html_path.resolve()),
                svg_path=str(svg_path.resolve()) if svg_path is not None else "",
                pdf_path=str(pdf_path.resolve()) if pdf_path is not None else "",
                status="ok",
            )
        )
        generated += 1
        outputs = [html_path.name]
        if pdf_path is not None:
            outputs.append(pdf_path.name)
        print("Generated:", contig, "->", ", ".join(outputs))

    for item in capped:
        manifest_rows.append(
            {
                "contig_id": item["contig_id"],
                "length": "",
                "mean_depth": "",
                "breadth_1x": "",
                "breadth_5x": "",
                "evidence_grade": (item["row"].get("evidence_grade") or "").strip(),
                "html_path": "",
                "svg_path": "",
                "pdf_path": "",
                "status": "skipped_cap",
            }
        )
    if capped:
        print("Capped: {} contigs beyond --max-contigs {} were not plotted".format(len(capped), args.max_contigs))

    manifest_path = outdir / "coverage_graphic_summary_manifest.tsv"
    write_manifest(manifest_path, manifest_rows)
    print("Manifest:", manifest_path.resolve())


if __name__ == "__main__":
    main()
