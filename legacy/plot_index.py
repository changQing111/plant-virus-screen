#!/usr/bin/env python3
"""Plot index page for a plant-virus-screen report directory.

Scans the ``nt_plots/`` and ``protein_plots/`` subdirectories (batch
manifests first, filename fallback second) and writes a self-contained
``index.html`` at the report root that lists every generated figure with
relative links to the HTML view and the PDF export, plus contig/HSP
counts, so figures can be browsed and located quickly.

Used by report.py after plot generation and runnable standalone:

    python3 plot_index.py --report-dir <report_dir>
"""

import argparse
import csv
import html
from pathlib import Path

NT_DIR = "nt_plots"
PROTEIN_DIR = "protein_plots"
NT_MANIFEST = "reference_graphic_summary_manifest.tsv"
PROTEIN_MANIFEST = "protein_graphic_summary_manifest.tsv"


def _read_manifest(path):
    rows = []
    if not path.is_file():
        return rows
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(row)
    return rows


def _scan_fallback(plots_dir):
    """Manifest-free fallback: derive entries from figure filenames."""
    entries = []
    for path in sorted(plots_dir.glob("*.graphic_summary.html")):
        stem = path.name[: -len(".graphic_summary.html")]
        name = stem.replace("__", " ").replace("_", " ")
        pdf_path = path.with_suffix(".pdf")
        entries.append(
            {
                "name": name,
                "contigs": "",
                "hsps": "",
                "html": path.name,
                "pdf": pdf_path.name if pdf_path.is_file() else "",
                "status": "ok",
            }
        )
    return entries


def collect_entries(report_dir, subdir, manifest_name, key_field, label_field):
    plots_dir = Path(report_dir) / subdir
    if not plots_dir.is_dir():
        return None
    rows = _read_manifest(plots_dir / manifest_name)
    if rows:
        entries = []
        for row in rows:
            html_name = Path(row.get("html_path") or "").name
            pdf_field = (row.get("pdf_path") or "").strip()
            pdf_name = Path(pdf_field).name if pdf_field else ""
            if pdf_name and not (plots_dir / pdf_name).is_file():
                pdf_name = ""
            entries.append(
                {
                    "name": (row.get(label_field) or row.get(key_field) or "").strip(),
                    "contigs": (row.get("contig_count") or "").strip(),
                    "hsps": (row.get("hsp_count") or "").strip(),
                    "html": html_name if html_name and (plots_dir / html_name).is_file() else "",
                    "pdf": pdf_name,
                    "status": (row.get("status") or "").strip(),
                }
            )
        return entries
    return _scan_fallback(plots_dir)


def _section_html(title, subdir, entries):
    parts = []
    generated = [e for e in entries if e["html"]]
    parts.append("<h2>{}（{}/，共 {} 张）</h2>".format(html.escape(title), html.escape(subdir), len(generated)))
    parts.append("<table>")
    parts.append("<tr><th>#</th><th>名称</th><th>contigs</th><th>HSPs</th><th>查看（HTML）</th><th>PDF</th></tr>")
    shown = 0
    for entry in entries:
        if not entry["html"]:
            if entry["status"]:
                parts.append(
                    "<tr class='muted'><td>-</td><td>{}</td><td colspan='4'>{}</td></tr>".format(
                        html.escape(entry["name"]), html.escape(entry["status"])
                    )
                )
            continue
        shown += 1
        pdf_cell = (
            "<a href='{}'>PDF</a>".format(html.escape(subdir + "/" + entry["pdf"])) if entry["pdf"] else ""
        )
        parts.append(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td><a href='{}'>HTML</a></td><td>{}</td></tr>".format(
                shown,
                html.escape(entry["name"]),
                html.escape(entry["contigs"]),
                html.escape(entry["hsps"]),
                html.escape(subdir + "/" + entry["html"]),
                pdf_cell,
            )
        )
    parts.append("</table>")
    return "\n".join(parts)


def build_index(report_dir, out_name="index.html"):
    """Write ``index.html`` into ``report_dir``; returns its path or None."""
    report_dir = Path(report_dir)
    section_specs = [
        ("nt 参考序列比对图", NT_DIR, NT_MANIFEST, "reference", "reference_label"),
        ("蛋白参考比对图", PROTEIN_DIR, PROTEIN_MANIFEST, "protein", "protein_label"),
    ]
    sections = []
    figure_count = 0
    for title, subdir, manifest, key_field, label_field in section_specs:
        entries = collect_entries(report_dir, subdir, manifest, key_field, label_field)
        if entries:
            sections.append(_section_html(title, subdir, entries))
            figure_count += sum(1 for entry in entries if entry["html"])
    if not sections:
        return None
    parts = []
    parts.append("<!DOCTYPE html>")
    parts.append("<html lang='zh-CN'>")
    parts.append("<head>")
    parts.append("<meta charset='utf-8'>")
    parts.append("<title>图形索引</title>")
    parts.append(
        "<style>"
        "body{font-family:Arial,'Microsoft YaHei',Helvetica,sans-serif;margin:24px;color:#222;background:#fff;}"
        "h1{font-size:22px;margin:0 0 8px 0;}"
        "h2{font-size:17px;margin:22px 0 10px 0;}"
        "p{font-size:14px;color:#444;}"
        "table{border-collapse:collapse;font-size:13px;margin:8px 0;}"
        "th,td{border:1px solid #d0d0d0;padding:4px 10px;text-align:left;}"
        "th{background:#f0f0f0;}"
        "a{color:#185fa5;text-decoration:none;}"
        "a:hover{text-decoration:underline;}"
        ".muted{color:#888;background:#fafafa;}"
        "</style>"
    )
    parts.append("</head>")
    parts.append("<body>")
    parts.append("<h1>图形索引</h1>")
    parts.append("<p>共 {} 张图形摘要，链接均相对本页。HTML 用于浏览器查看（含悬停提示），PDF 为矢量导出，可缩放不失真。</p>".format(figure_count))
    if (report_dir / "REPORT.md").is_file():
        parts.append("<p><a href='REPORT.md'>返回报告（REPORT.md）</a></p>")
    parts.extend(sections)
    parts.append("</body>")
    parts.append("</html>")
    out_path = report_dir / out_name
    out_path.write_text("\n".join(parts), encoding="utf-8")
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate the figure index page (index.html) for a report directory."
    )
    parser.add_argument(
        "--report-dir",
        required=True,
        help="Report directory containing nt_plots/ and/or protein_plots/",
    )
    args = parser.parse_args()
    out_path = build_index(args.report_dir)
    if out_path is None:
        raise SystemExit("No plot directories found under {}; nothing to index.".format(args.report_dir))
    print("Wrote index:", out_path.resolve())


if __name__ == "__main__":
    main()
