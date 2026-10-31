#!/usr/bin/env python3

import argparse
import csv
from collections import OrderedDict
from pathlib import Path


def split_sorted_join(values):
    cleaned = sorted(v for v in values if v)
    return ";".join(cleaned)


def load_rows(path):
    with Path(path).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize_rows(rows):
    grouped = OrderedDict()

    for row in rows:
        virus = (row.get("putative_virus") or "").strip()
        if not virus:
            virus = "Unknown virus"

        if virus not in grouped:
            grouped[virus] = {
                "putative_virus": virus,
                "genome_composition": row.get("genome_composition", "").strip(),
                "virus_genome_type_cn": row.get("virus_genome_type_cn", "").strip(),
                "baltimore_group": row.get("baltimore_group", "").strip(),
                "host_source": row.get("host_source", "").strip(),
                "plant_virus_status": row.get("plant_virus_status", "").strip(),
                "contig_ids": set(),
                "group_ids": set(),
                "candidate_names": set(),
                "segment_or_protein": set(),
                "priority_counts": {"primary": 0, "secondary": 0, "other": 0},
            }

        bucket = grouped[virus]

        rep_id = (row.get("representative_query_id") or "").strip()
        if rep_id:
            bucket["contig_ids"].add(rep_id)

        group_id = (row.get("group_id") or "").strip()
        if group_id:
            bucket["group_ids"].add(group_id)

        candidate_name = (row.get("candidate_name") or "").strip()
        if candidate_name:
            bucket["candidate_names"].add(candidate_name)

        segment = (row.get("putative_segment_or_protein") or "").strip()
        if segment:
            bucket["segment_or_protein"].add(segment)

        priority = (row.get("priority") or "").strip().lower()
        if priority in bucket["priority_counts"]:
            bucket["priority_counts"][priority] += 1
        else:
            bucket["priority_counts"]["other"] += 1

        if not bucket["genome_composition"]:
            bucket["genome_composition"] = (row.get("genome_composition") or "").strip()
        if not bucket["virus_genome_type_cn"]:
            bucket["virus_genome_type_cn"] = (row.get("virus_genome_type_cn") or "").strip()
        if not bucket["baltimore_group"]:
            bucket["baltimore_group"] = (row.get("baltimore_group") or "").strip()
        if not bucket["host_source"]:
            bucket["host_source"] = (row.get("host_source") or "").strip()
        if not bucket["plant_virus_status"]:
            bucket["plant_virus_status"] = (row.get("plant_virus_status") or "").strip()

    summary = []
    for bucket in grouped.values():
        summary.append(
            {
                "putative_virus": bucket["putative_virus"],
                "genome_composition": bucket["genome_composition"],
                "virus_genome_type_cn": bucket["virus_genome_type_cn"],
                "baltimore_group": bucket["baltimore_group"],
                "host_source": bucket["host_source"],
                "plant_virus_status": bucket["plant_virus_status"] or "unknown",
                "supporting_contig_count": str(len(bucket["contig_ids"])),
                "supporting_group_count": str(len(bucket["group_ids"])),
                "primary_count": str(bucket["priority_counts"]["primary"]),
                "secondary_count": str(bucket["priority_counts"]["secondary"]),
                "segment_or_protein": split_sorted_join(bucket["segment_or_protein"]),
                "candidate_names": split_sorted_join(bucket["candidate_names"]),
                "contig_ids": split_sorted_join(bucket["contig_ids"]),
            }
        )

    summary.sort(
        key=lambda row: (
            -int(row["supporting_contig_count"]),
            -int(row["primary_count"]),
            row["putative_virus"].lower(),
        )
    )
    return summary


def write_summary(path, summary_rows):
    fieldnames = [
        "putative_virus",
        "genome_composition",
        "virus_genome_type_cn",
        "baltimore_group",
        "host_source",
        "plant_virus_status",
        "supporting_contig_count",
        "supporting_group_count",
        "primary_count",
        "secondary_count",
        "segment_or_protein",
        "candidate_names",
        "contig_ids",
    ]

    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize final_retained_candidates.tsv by virus species and supporting contig counts."
    )
    parser.add_argument(
        "input_tsv",
        nargs="?",
        default="final_retained_candidates.tsv",
        help="Input final_retained_candidates.tsv file",
    )
    parser.add_argument(
        "-o",
        "--out",
        default="virus_species_summary.tsv",
        help="Output summary TSV file",
    )
    args = parser.parse_args()

    rows = load_rows(args.input_tsv)
    summary_rows = summarize_rows(rows)
    write_summary(args.out, summary_rows)


if __name__ == "__main__":
    main()
