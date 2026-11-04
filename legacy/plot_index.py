#!/usr/bin/env python3
"""Plot index page for a plant-virus-screen report directory.

Scans the ``nt_plots/``, ``protein_plots/`` and ``coverage_plots``
subdirectories (batch manifests first, filename fallback second) and
writes a self-contained ``index.html`` at the report root that lists
every generated figure with relative links to the HTML view and the PDF
export, plus per-kind counts (contigs/HSPs for alignment figures,
length/mean depth/breadth for coverage figures), so figures can be
browsed and located quickly.

Used by report.py after plot generation and runnable standalone:

    python3 plot_index.py --report-dir <report_dir>
"""

import argparse
import csv
import html
from pathlib import Path

NT_DIR = "nt_plots"
PROTEIN_DIR = "protein_plots"
COVERAGE_DIR = "coverage_plots"
NT_MANIFEST = "reference_graphic_summary_manifest.tsv"
PROTEIN_MANIFEST = "protein_graphic_summary_manifest.tsv"
COVERAGE_MANIFEST = "coverage_graphic_summary_manifest.tsv"


def _read_manifest(path):
    rows = []
    if not path.is_file():
        return rows
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(row)
    return rows


def _scan_fallback(plots_dir, pattern, suffix):
    """Manifest-free fallback: derive entries from figure filenames."""
    entries = []
    for path in sorted(plots_dir.glob(pattern)):
        stem = path.name[: -len(suffix)]
        name = stem.replace("__", " ").replace("_", " ")
        pdf_path = path.with_suffix(".pdf")
        entries.append(
            {
                "name": name,
                "contigs": "",
                "hsps": "",
                "length": "",
                "mean_depth": "",
                "breadth_1x": "",
                "html": path.name,
                "pdf": pdf_path.name if pdf_path.is_file() else "",
                "status": "ok",
            }
        )
    return entries


def collect_entries(report_dir, spec):
    plots_dir = Path(report_dir) / spec["subdir"]
    if not plots_dir.is_dir():
        return None
    rows = _read_manifest(plots_dir / spec["manifest"])
    if rows:
        entries = []
        for row in rows:
            html_name = Path(row.get("html_path") or "").name
            pdf_field = (row.get("pdf_path") or "").strip()
            pdf_name = Path(pdf_field).name if pdf_field else ""
            if pdf_name and not (plots_dir / pdf_name).is_file():
                pdf_name = ""
            entry = {
                "name": (row.get(spec["label_field"]) or row.get(spec["key_field"]) or "").strip(),
                "contigs": (row.get("contig_count") or "").strip(),
                "hsps": (row.get("hsp_count") or "").strip(),
                "length": (row.get("length") or "").strip(),
                "mean_depth": (row.get("mean_depth") or "").strip(),
                "breadth_1x": (row.get("breadth_1x") or "").strip(),
                "html": html_name if html_name and (plots_dir / html_name).is_file() else "",
                "pdf": pdf_name,
                "status": (row.get("status") or "").strip(),
            }
            entries.append(entry)
        return entries
    return _scan_fallback(plots_dir, spec["glob"], spec["suffix"])


def _section_html(title, subdir, entries, kind):
    parts = []
    generated = [e for e in entries if e["html"]]
    parts.append("<h2>{}（{}/，共 {} 张）</h2>".format(html.escape(title), html.escape(subdir), len(generated)))
    parts.append("<table>")
    if kind == "coverage":
        header = "<tr><th>#</th><th>contig</th><th>长度 (bp)</th><th>平均深度</th><th>≥1x 广度</th><th>查看（HTML）</th><th>PDF</th></tr>"
        span = 5
    else:
        header = "<tr><th>#</th><th>名称</th><th>contigs</th><th>HSPs</th><th>查看（HTML）</th><th>PDF</th></tr>"
        span = 4
    parts.append(header)
    shown = 0
    for entry in entries:
        if not entry["html"]:
            if entry["status"]:
                parts.append(
                    "<tr class='muted'><td>-</td><td>{}</td><td colspan='{}'>{}</td></tr>".format(
                        html.escape(entry["name"]), span, html.escape(entry["status"])
                    )
                )
            continue
        shown += 1
        pdf_cell = (
            "<a href='{}'>PDF</a>".format(html.escape(subdir + "/" + entry["pdf"])) if entry["pdf"] else ""
        )
        link_cell = "<a href='{}'>HTML</a>".format(html.escape(subdir + "/" + entry["html"]))
        if kind == "coverage":
            parts.append(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
                    shown,
                    html.escape(entry["name"]),
                    html.escape(entry["length"]),
                    html.escape(entry["mean_depth"]),
                    html.escape(entry["breadth_1x"]),
                    link_cell,
                    pdf_cell,
                )
            )
        else:
            parts.append(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
                    shown,
                    html.escape(entry["name"]),
                    html.escape(entry["contigs"]),
                    html.escape(entry["hsps"]),
                    link_cell,
                    pdf_cell,
                )
            )
    parts.append("</table>")
    return "\n".join(parts)


def _section_specs():
    return [
        {
            "title": "nt 参考序列比对图",
            "subdir": NT_DIR,
            "manifest": NT_MANIFEST,
            "key_field": "reference",
            "label_field": "reference_label",
            "kind": "alignment",
            "glob": "*.graphic_summary.html",
            "suffix": ".graphic_summary.html",
        },
        {
            "title": "蛋白参考比对图",
            "subdir": PROTEIN_DIR,
            "manifest": PROTEIN_MANIFEST,
            "key_field": "protein",
            "label_field": "protein_label",
            "kind": "alignment",
            "glob": "*.graphic_summary.html",
            "suffix": ".graphic_summary.html",
        },
        {
            "title": "候选覆盖深度图",
            "subdir": COVERAGE_DIR,
            "manifest": COVERAGE_MANIFEST,
            "key_field": "contig_id",
            "label_field": "contig_id",
            "kind": "coverage",
            "glob": "*.coverage_summary.html",
            "suffix": ".coverage_summary.html",
        },
    ]


def build_index(report_dir, out_name="index.html"):
    """Write ``index.html`` into ``report_dir``; returns its path or None."""
    report_dir = Path(report_dir)
    sections = []
    figure_count = 0
    for spec in _section_specs():
        entries = collect_entries(report_dir, spec)
        if entries:
            sections.append(_section_html(spec["title"], spec["subdir"], entries, spec["kind"]))
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
        help="Report directory containing nt_plots/, protein_plots/ and/or coverage_plots/",
    )
    args = parser.parse_args()
    out_path = build_index(args.report_dir)
    if out_path is None:
        raise SystemExit("No plot directories found under {}; nothing to index.".format(args.report_dir))
    print("Wrote index:", out_path.resolve())


if __name__ == "__main__":
    main()
