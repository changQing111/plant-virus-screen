#!/usr/bin/env python3

import argparse
import csv
import importlib.util
import re
from pathlib import Path


def load_visualizer_module(module_path):
    spec = importlib.util.spec_from_file_location("visualize_blastx_protein", str(module_path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch-generate blastx protein graphic summaries (HTML + SVG + PDF) from final_retained_candidates.tsv."
    )
    parser.add_argument(
        "--final-tsv",
        default="final_retained_candidates.tsv",
        help="Input final_retained_candidates.tsv",
    )
    parser.add_argument(
        "--blastx",
        default="diamond_blastx_res.txt",
        help="Input diamond blastx result file",
    )
    parser.add_argument(
        "--outdir",
        default="protein_graphic_summaries",
        help="Output directory for generated HTML/PDF files",
    )
    parser.add_argument(
        "--module",
        default="visualize_blastx_protein.py",
        help="Path to visualize_blastx_protein.py",
    )
    parser.add_argument(
        "--protein-fasta",
        default="/home/changq/my_data/Document/viral/uniprotkb_taxonomy_id_10239_2026_04_03.fasta",
        help="Protein FASTA used to resolve protein labels.",
    )
    parser.add_argument(
        "--top-contigs",
        type=int,
        default=0,
        help="Optional cap on the number of contigs displayed in each figure (0 means all).",
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
        help="Skip PDF export.",
    )
    parser.add_argument(
        "--no-svg",
        action="store_true",
        help="Skip separate SVG export.",
    )
    return parser.parse_args()


def extract_proteins(final_tsv):
    proteins = []
    seen = set()
    with Path(final_tsv).open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            best_protein_hit = (row.get("best_protein_hit") or "").strip()
            if not best_protein_hit or best_protein_hit.startswith("no strong protein hit"):
                continue
            protein = best_protein_hit.split(" | ", 1)[0].strip()
            if protein and protein not in seen:
                seen.add(protein)
                proteins.append(protein)
    return proteins


def protein_filename_prefix(protein):
    parts = protein.split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1]
    return protein


def protein_file_label(protein_label, protein):
    if "OS=" in protein_label:
        name = protein_label.split("OS=", 1)[1]
        name = re.split(r"\s+OX=|\s+GN=|\s+PE=|\s+SV=|$", name, maxsplit=1)[0]
    elif " | " in protein_label:
        name = protein_label.split(" | ", 1)[1]
    else:
        name = protein_label
    slug = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
    prefix = protein_filename_prefix(protein)
    if not slug:
        slug = prefix
    return "{}__{}".format(prefix, slug)


def limit_grouped_hits(module, grouped_hits, top_contigs):
    if top_contigs <= 0:
        return grouped_hits
    sorted_names = sorted(
        grouped_hits,
        key=lambda name: (
            module.anchor_hit(grouped_hits[name]).subject_min,
            module.best_bitscore(grouped_hits[name]),
            max(hit.qlen for hit in grouped_hits[name]),
        ),
    )
    sorted_names = sorted_names[:top_contigs]
    return {name: grouped_hits[name] for name in sorted_names}


def write_manifest(path, rows):
    fieldnames = [
        "protein",
        "protein_label",
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
    proteins = extract_proteins(args.final_tsv)
    manifest_rows = []

    for protein in proteins:
        hits = module.load_blastx_hits(args.blastx, protein, args.min_align_aa)
        allowed = module.load_allowed_contigs(args.final_tsv, protein)
        if allowed is not None:
            hits = [hit for hit in hits if hit.qseqid in allowed]

        if not hits:
            manifest_rows.append(
                {
                    "protein": protein,
                    "protein_label": protein,
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
        protein_label = module.load_protein_label(protein, args.protein_fasta, args.final_tsv)

        prefix = module.sanitize_filename(protein_file_label(protein_label, protein))
        html_path = outdir / (prefix + ".graphic_summary.html")
        html_text = module.build_html(protein, protein_label, grouped_hits)
        html_path, svg_path, pdf_path = module.write_output_bundle(
            html_path,
            html_text,
            write_svg=not args.no_svg,
            write_pdf=not args.no_pdf,
        )

        manifest_rows.append(
            {
                "protein": protein,
                "protein_label": protein_label,
                "html_path": str(html_path.resolve()),
                "svg_path": str(svg_path.resolve()) if svg_path is not None else "",
                "pdf_path": str(pdf_path.resolve()) if pdf_path is not None else "",
                "contig_count": str(len(grouped_hits)),
                "hsp_count": str(sum(len(hs) for hs in grouped_hits.values())),
                "status": "ok",
            }
        )

        outputs = [html_path.name]
        if svg_path is not None:
            outputs.append(svg_path.name)
        if pdf_path is not None:
            outputs.append(pdf_path.name)
        print("Generated:", protein, "->", ", ".join(outputs))

    manifest_path = outdir / "protein_graphic_summary_manifest.tsv"
    write_manifest(manifest_path, manifest_rows)
    print("Manifest:", manifest_path.resolve())


if __name__ == "__main__":
    main()
