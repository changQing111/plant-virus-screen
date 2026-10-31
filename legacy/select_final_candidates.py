#!/usr/bin/env python3

import argparse
import csv
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


ORTHOTOSPO_KEYWORDS = (
    "orthotospovirus",
    "tomato zonate spot virus",
    "calla lily chlorotic spot virus",
    "tomato necrotic spot",
    "tospovirus kiwifruit",
    "mulberry vein banding virus",
    "alstroemeria yellow spot virus",
    "watermelon silver mottle virus",
)

SECONDARY_EXCLUDE_KEYWORDS = (
    "partitivirus",
    "ourmia-like",
    "caulimovirus",
    "retro-like virus",
    "retrovirus",
    "megavirus",
    "mimivirus",
    "pandoravirus",
    "marseillevirus",
    "klosneuvirus",
    "tupanvirus",
    "medusavirus",
    "fadolivirus",
    "poxvirus",
    "crovi",
    "crovb",
    "phage",
    "caudovir",
    "siphovir",
    "myovir",
    "podovir",
)

LABEL_PRIORITY = {
    "dual_support_divergent": 4,
    "dual_support_close": 3,
    "nt_only_strong": 2,
    "aa_only_strong": 2,
    "suspect_protein_only": 1,
    "no_strong_hit": 0,
}

XLSX_NS = {
    "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


def safe_float(value):
    if value is None:
        return 0.0
    text = str(value).strip()
    if not text:
        return 0.0
    return float(text)


def load_rows(path):
    rows = []
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def annotation_text(row):
    text = " ".join(
        [
            row.get("blastn_hit", ""),
            row.get("blastn_anno", ""),
            row.get("blastx_hit", ""),
            row.get("blastx_anno", ""),
            row.get("suspect_blastx_hit", ""),
            row.get("suspect_blastx_anno", ""),
        ]
    )
    return text.lower()


def has_any(text, keywords):
    return any(keyword in text for keyword in keywords)


def is_orthotospovirus_like(row):
    return has_any(annotation_text(row), ORTHOTOSPO_KEYWORDS)


def role_l_rdrp(row):
    return "rna-directed rna polymerase l" in row.get("blastx_anno", "").lower()


def role_m_glycoprotein(row):
    return "envelopment polyprotein" in row.get("blastx_anno", "").lower()


def role_m_nsm(row):
    return re.search(r"\bnsm\b", row.get("blastx_anno", "").lower()) is not None


def role_s_nss(row):
    anno = row.get("blastx_anno", "").lower()
    return (
        "non-structural protein" in anno or "nonstructural protein" in anno
    ) and "nucleoprotein" not in anno


def role_s_n(row):
    return "nucleoprotein" in row.get("blastx_anno", "").lower()


def rank_row_supported(row):
    return (
        LABEL_PRIORITY.get(row.get("label", ""), 0),
        safe_float(row.get("combined_score")),
        safe_float(row.get("blastx_qcovhsp")),
        safe_float(row.get("blastn_qcovs")),
        safe_float(row.get("blastx_align_len_aa")),
        safe_float(row.get("blastn_align_len")),
        safe_float(row.get("blastx_pident")),
        safe_float(row.get("blastn_pident")),
        safe_float(row.get("qlen")),
        safe_float(row.get("qcov_from_id")),
    )


def rank_row_longest(row):
    return (
        LABEL_PRIORITY.get(row.get("label", ""), 0),
        safe_float(row.get("combined_score")),
        safe_float(row.get("qlen")),
        safe_float(row.get("blastx_align_len_aa")),
        safe_float(row.get("blastn_align_len")),
        safe_float(row.get("blastx_qcovhsp")),
        safe_float(row.get("blastn_qcovs")),
        safe_float(row.get("qcov_from_id")),
    )


def rank_row_target_3000(row):
    qlen = safe_float(row.get("qlen"))
    return (
        LABEL_PRIORITY.get(row.get("label", ""), 0),
        safe_float(row.get("combined_score")),
        -abs(qlen - 3000.0),
        safe_float(row.get("blastx_align_len_aa")),
        safe_float(row.get("blastn_align_len")),
        safe_float(row.get("blastx_qcovhsp")),
        safe_float(row.get("blastn_qcovs")),
        qlen,
        safe_float(row.get("qcov_from_id")),
    )


def choose_best(rows, predicate, rank_fn):
    selected = [row for row in rows if predicate(row)]
    if not selected:
        return None
    selected.sort(key=rank_fn, reverse=True)
    return selected[0]


def choose_top_n(rows, predicate, rank_fn, n, exclude_ids=None):
    exclude_ids = exclude_ids or set()
    selected = [row for row in rows if predicate(row) and row.get("query_id") not in exclude_ids]
    if not selected:
        return []
    selected.sort(key=rank_fn, reverse=True)
    return selected[:n]


def cell_ref_to_index(cell_ref):
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    index = 0
    for ch in letters:
        index = index * 26 + (ord(ch.upper()) - ord("A") + 1)
    return index - 1


def xlsx_cell_value(cell, shared_strings):
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(t.text or "" for t in cell.iterfind(".//a:t", XLSX_NS))
    value_node = cell.find("a:v", XLSX_NS)
    if value_node is None:
        return ""
    value = value_node.text or ""
    if cell_type == "s" and value:
        return shared_strings[int(value)]
    return value


def iter_xlsx_rows(path, sheet_name):
    with zipfile.ZipFile(path) as archive:
        shared_strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for si in shared_root.findall("a:si", XLSX_NS):
                shared_strings.append("".join(t.text or "" for t in si.iterfind(".//a:t", XLSX_NS)))

        workbook_root = ET.fromstring(archive.read("xl/workbook.xml"))
        rels_root = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels_root}

        target = None
        for sheet in workbook_root.find("a:sheets", XLSX_NS):
            if sheet.attrib.get("name") == sheet_name:
                rid = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
                target = "xl/" + rel_map[rid]
                break
        if target is None:
            raise ValueError("Sheet not found: %s" % sheet_name)

        sheet_root = ET.fromstring(archive.read(target))
        for row in sheet_root.findall(".//a:sheetData/a:row", XLSX_NS):
            values = {}
            max_index = -1
            for cell in row.findall("a:c", XLSX_NS):
                ref = cell.attrib.get("r", "")
                index = cell_ref_to_index(ref)
                values[index] = xlsx_cell_value(cell, shared_strings)
                max_index = max(max_index, index)
            if max_index < 0:
                yield []
                continue
            yield [values.get(i, "") for i in range(max_index + 1)]


def normalize_name(text):
    text = (text or "").lower().strip()
    text = text.replace("_", " ")
    # Fold punctuation variants such as "spot-associated" and "spot associated"
    # into the same token stream before matching against VMR names.
    text = re.sub(r"[-/()\\[\\],;:+]+", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def split_multi_value(text):
    text = (text or "").strip()
    if not text:
        return []
    parts = re.split(r"\s*;\s*|\s*,\s*", text)
    return [part.strip() for part in parts if part.strip()]


def extract_accessions(text):
    text = text or ""
    accessions = set()
    for token in re.findall(r"[A-Z]{1,4}_?\d+(?:\.\d+)?", text):
        accessions.add(token)
        accessions.add(token.split(".")[0])
    return list(accessions)


def extract_os_name(text):
    text = text or ""
    match = re.search(r"\bOS=(.+?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)", text)
    if match:
        return match.group(1).strip()
    return ""


def genome_type_cn(genome):
    mapping = {
        "dsdna": "双链DNA病毒",
        "ssdna": "单链DNA病毒",
        "dsrna": "双链RNA病毒",
        "ssrna(+)": "正单链RNA病毒",
        "ssrna(-)": "负单链RNA病毒",
        "ssrna(+/-)": "负单链RNA病毒（双义编码）",
        "ssrna(-/+)": "负单链RNA病毒（双义编码）",
        "ssrna": "单链RNA病毒",
        "ssrna-rt": "单链RNA逆转录病毒",
        "dsdna-rt": "双链DNA逆转录病毒",
    }
    return mapping.get((genome or "").lower(), "")


def baltimore_group(genome):
    mapping = {
        "dsdna": "I",
        "ssdna": "II",
        "dsrna": "III",
        "ssrna(+)": "IV",
        "ssrna(-)": "V",
        "ssrna(+/-)": "V",
        "ssrna(-/+)": "V",
        "ssrna-rt": "VI",
        "dsdna-rt": "VII",
    }
    return mapping.get((genome or "").lower(), "")


def add_consistent(mapping, key, genome):
    if not key or not genome:
        return
    mapping.setdefault(key, set()).add(genome)


def unique_mapping_value(mapping, key):
    values = mapping.get(key, set())
    if len(values) == 1:
        return list(values)[0]
    return ""


def load_vmr_genome_maps(path):
    rows = list(iter_xlsx_rows(path, "VMR MSL41"))
    if not rows:
        raise ValueError("No rows found in VMR sheet")
    header = rows[0]
    index = {name: i for i, name in enumerate(header)}

    accession_map = {}
    name_map = {}
    genus_map = {}
    species_map = {}
    accession_host_map = {}
    name_host_map = {}
    genus_host_map = {}
    species_host_map = {}

    for row in rows[1:]:
        if not row:
            continue
        genome = row[index["Genome"]].strip() if "Genome" in index and index["Genome"] < len(row) else ""

        virus_names = row[index["Virus name(s)"]] if "Virus name(s)" in index and index["Virus name(s)"] < len(row) else ""
        virus_abbrs = row[index["Virus name abbreviation(s)"]] if "Virus name abbreviation(s)" in index and index["Virus name abbreviation(s)"] < len(row) else ""
        genus = row[index["Genus"]] if "Genus" in index and index["Genus"] < len(row) else ""
        species = row[index["Species"]] if "Species" in index and index["Species"] < len(row) else ""
        accessions = row[index["Virus GENBANK accession"]] if "Virus GENBANK accession" in index and index["Virus GENBANK accession"] < len(row) else ""
        host_source = row[index["Host source"]].strip() if "Host source" in index and index["Host source"] < len(row) else ""

        if genome:
            for acc in extract_accessions(accessions):
                add_consistent(accession_map, acc, genome)
            for name in split_multi_value(virus_names):
                add_consistent(name_map, normalize_name(name), genome)
            for abbr in split_multi_value(virus_abbrs):
                add_consistent(name_map, normalize_name(abbr), genome)
            add_consistent(genus_map, normalize_name(genus), genome)
            add_consistent(species_map, normalize_name(species), genome)

        if host_source:
            for acc in extract_accessions(accessions):
                add_consistent(accession_host_map, acc, host_source)
            for name in split_multi_value(virus_names):
                add_consistent(name_host_map, normalize_name(name), host_source)
            for abbr in split_multi_value(virus_abbrs):
                add_consistent(name_host_map, normalize_name(abbr), host_source)
            add_consistent(genus_host_map, normalize_name(genus), host_source)
            add_consistent(species_host_map, normalize_name(species), host_source)

    return (
        accession_map,
        name_map,
        genus_map,
        species_map,
        accession_host_map,
        name_host_map,
        genus_host_map,
        species_host_map,
    )


def infer_genome_from_putative_virus(putative_virus, genus_map):
    value = (putative_virus or "").strip()
    if not value:
        return ""
    key = normalize_name(value.replace("-like", "").replace("_like", ""))
    return unique_mapping_value(genus_map, key)


def classify_plant_virus(host_source):
    text = (host_source or "").strip().lower()
    if not text:
        return "unknown"
    plant_keywords = (
        "plant",
        "plants",
        "viridiplantae",
        "green plant",
        "angiosperm",
        "monocot",
        "dicot",
    )
    if any(keyword in text for keyword in plant_keywords):
        return "plant virus"
    return "non-plant virus"


def annotate_with_vmr(rows, vmr_xlsx):
    (
        accession_map,
        name_map,
        genus_map,
        species_map,
        accession_host_map,
        name_host_map,
        genus_host_map,
        species_host_map,
    ) = load_vmr_genome_maps(vmr_xlsx)
    for row in rows:
        genome = ""
        host_source = ""

        candidate_strings = [
            row.get("best_nt_hit", ""),
            row.get("best_protein_hit", ""),
            row.get("putative_virus", ""),
            row.get("candidate_name", ""),
        ]

        for text in candidate_strings:
            if genome:
                break
            for acc in extract_accessions(text):
                genome = unique_mapping_value(accession_map, acc)
                if not host_source:
                    host_source = unique_mapping_value(accession_host_map, acc)
                if genome:
                    break

        names_to_try = []
        nt_hit = row.get("best_nt_hit", "")
        if " | " in nt_hit:
            names_to_try.append(nt_hit.split(" | ", 1)[1].split(", complete sequence")[0].strip())
        protein_hit = row.get("best_protein_hit", "")
        if " | " in protein_hit:
            names_to_try.append(protein_hit.split(" | ", 1)[1])
        os_name = extract_os_name(protein_hit)
        if os_name:
            names_to_try.append(os_name)
        names_to_try.append(row.get("putative_virus", "").replace("-like", ""))

        for name in names_to_try:
            norm = normalize_name(name)
            if not genome:
                genome = (
                    unique_mapping_value(name_map, norm)
                    or unique_mapping_value(species_map, norm)
                    or unique_mapping_value(genus_map, norm)
                )
            if not host_source:
                host_source = (
                    unique_mapping_value(name_host_map, norm)
                    or unique_mapping_value(species_host_map, norm)
                    or unique_mapping_value(genus_host_map, norm)
                )
            if genome and host_source:
                break

        if not genome:
            genome = infer_genome_from_putative_virus(row.get("putative_virus", ""), genus_map)
        if not host_source:
            host_source = unique_mapping_value(
                genus_host_map,
                normalize_name(row.get("putative_virus", "").replace("-like", "").replace("_like", "")),
            )

        row["genome_composition"] = genome
        row["virus_genome_type_cn"] = genome_type_cn(genome)
        row["baltimore_group"] = baltimore_group(genome)
        row["host_source"] = host_source
        row["plant_virus_status"] = classify_plant_virus(host_source)

    return rows


def format_hit(hit_id, hit_anno, empty_message):
    hit_id = (hit_id or "").strip()
    hit_anno = (hit_anno or "").strip()
    if not hit_id:
        return empty_message
    if hit_anno:
        return "%s | %s" % (hit_id, hit_anno)
    return hit_id


def family_from_secondary(row):
    text = annotation_text(row)
    if "potyvirus" in text or "plum pox virus" in text or "_9potv" in text:
        return "Potyvirus-like"
    if "ourmia-like" in text:
        return "Ourmia-like"
    if "partitivirus" in text:
        return "Partitivirus-like"
    if "caulimovirus" in text:
        return "Caulimovirus-like"
    return "Other_virus_like"


def infer_putative_virus_name(row):
    protein_hit = row.get("best_protein_hit", "")
    os_name = extract_os_name(protein_hit)
    if os_name:
        return os_name

    blastx_anno = row.get("blastx_anno", "")
    os_name = extract_os_name(blastx_anno)
    if os_name:
        return os_name

    blastn_anno = row.get("blastn_anno", "")
    if blastn_anno:
        name = blastn_anno
        name = re.split(r"\s+(?:isolate|strain)\s+", name, maxsplit=1)[0]
        name = re.split(r"\s+segment\s+[A-Za-z0-9-]+", name, maxsplit=1)[0]
        name = re.split(r",\s*(?:complete|coding-complete|partial)\b", name, maxsplit=1)[0]
        return name.strip(" ,;")

    family = family_from_secondary(row)
    if family != "Other_virus_like":
        return family
    return "Unclassified virus-like candidate"


def infer_segment_or_protein(row):
    anno = row.get("blastx_anno", "").lower()
    if "rna-directed rna polymerase" in anno or "replicase" in anno:
        return "RdRP/replicase"
    if "nucleoprotein" in anno:
        return "Nucleoprotein"
    if "non-structural protein" in anno or "nonstructural protein" in anno:
        return "Non-structural protein"
    if "glycoprotein" in anno or "envelopment polyprotein" in anno:
        return "Glycoprotein"
    if "movement protein" in anno or re.search(r"\bnsm\b", anno):
        return "Movement protein"
    if "coat protein" in anno:
        return "Coat protein"
    if "polyprotein" in anno:
        return "Polyprotein"
    if row.get("blastn_hit"):
        return "Best nucleotide-supported candidate"
    return "Best protein-supported candidate"


def virus_group_key(row):
    key = normalize_name(infer_putative_virus_name(row))
    return key or normalize_name(row.get("query_id", ""))


def slugify(text, default_value):
    slug = re.sub(r"[^A-Za-z0-9]+", "_", text or "").strip("_")
    return slug or default_value


def primary_candidate_name(row):
    return slugify(infer_putative_virus_name(row), "Virus_candidate") + "_primary"


def backup_candidate_name(row, index):
    return slugify(infer_putative_virus_name(row), "Virus_candidate") + "_backup" + str(index)


def extra_secondary_candidate_name(row, index):
    return slugify(infer_putative_virus_name(row), "Virus_candidate") + "_secondary" + str(index)


def is_generic_primary_candidate(row):
    label = row.get("label", "")
    qlen = safe_float(row.get("qlen"))
    qcov_nt = safe_float(row.get("blastn_qcovs"))
    qcov_aa = safe_float(row.get("blastx_qcovhsp"))
    combined = safe_float(row.get("combined_score"))

    if label in ("no_strong_hit", "suspect_protein_only"):
        return False
    if label in ("dual_support_divergent", "dual_support_close"):
        return combined >= 300
    if label == "aa_only_strong":
        return qlen >= 1500 and qcov_aa >= 25
    if label == "nt_only_strong":
        return qlen >= 800 and qcov_nt >= 20
    return False


def is_generic_secondary_candidate(row, min_qlen, min_qcovhsp):
    label = row.get("label", "")
    qlen = safe_float(row.get("qlen"))
    qcov_nt = safe_float(row.get("blastn_qcovs"))
    qcov_aa = safe_float(row.get("blastx_qcovhsp"))

    if label in ("no_strong_hit", "suspect_protein_only"):
        return False
    if label in ("dual_support_divergent", "dual_support_close"):
        return qlen >= max(1000, min_qlen / 2.0)
    if label == "aa_only_strong":
        return qlen >= min_qlen and qcov_aa >= min_qcovhsp
    if label == "nt_only_strong":
        return qlen >= min_qlen and qcov_nt >= 20
    return False


def group_rows_by_virus(rows, predicate):
    grouped = {}
    for row in rows:
        if not predicate(row):
            continue
        grouped.setdefault(virus_group_key(row), []).append(row)
    return grouped


def generic_candidate_name(row, index):
    putative = infer_putative_virus_name(row)
    slug = re.sub(r"[^A-Za-z0-9]+", "_", putative).strip("_")
    if not slug:
        slug = "Generic_candidate"
    return "%s_%d" % (slug[:50], index)


def secondary_name(row):
    family = family_from_secondary(row).replace("-", "_")
    if "polyprotein" in row.get("blastx_anno", "").lower():
        return family + "_polyprotein"
    return family + "_candidate"


def secondary_segment(row):
    if "polyprotein" in row.get("blastx_anno", "").lower():
        return "Polyprotein"
    return "Secondary candidate"


def is_secondary_candidate(row, min_qlen, min_qcovhsp):
    text = annotation_text(row)
    if is_orthotospovirus_like(row):
        return False
    if row.get("label") != "aa_only_strong":
        return False
    if row.get("blastn_hit"):
        return False
    if has_any(text, SECONDARY_EXCLUDE_KEYWORDS):
        return False
    if "polyprotein" not in row.get("blastx_anno", "").lower():
        return False
    if safe_float(row.get("qlen")) < min_qlen:
        return False
    if safe_float(row.get("blastx_qcovhsp")) < min_qcovhsp:
        return False
    return True


def is_generic_fallback_candidate(row):
    text = annotation_text(row)
    label = row.get("label", "")
    qlen = safe_float(row.get("qlen"))
    qcov_nt = safe_float(row.get("blastn_qcovs"))
    qcov_aa = safe_float(row.get("blastx_qcovhsp"))
    combined = safe_float(row.get("combined_score"))

    if has_any(text, SECONDARY_EXCLUDE_KEYWORDS):
        return False
    if label in ("suspect_protein_only", "no_strong_hit"):
        return False
    if label == "dual_support_divergent":
        return combined >= 500 and qlen >= 1000
    if label == "aa_only_strong":
        return qlen >= 1500 and qcov_aa >= 25
    if label == "nt_only_strong":
        return qlen >= 800 and qcov_nt >= 20
    return False


def build_output_row(priority, candidate_name, putative_virus, segment, row, reason):
    return {
        "priority": priority,
        "candidate_name": candidate_name,
        "representative_query_id": row.get("query_id", ""),
        "group_id": row.get("group_id", ""),
        "putative_virus": putative_virus,
        "putative_segment_or_protein": segment,
        "status": "keep",
        "best_nt_hit": format_hit(
            row.get("blastn_hit", ""),
            row.get("blastn_anno", ""),
            "no strong nucleotide hit",
        ),
        "blastn_pident": row.get("blastn_pident", ""),
        "blastn_qcovs": row.get("blastn_qcovs", ""),
        "best_protein_hit": format_hit(
            row.get("blastx_hit", ""),
            row.get("blastx_anno", ""),
            "no strong protein hit",
        ),
        "blastx_pident": row.get("blastx_pident", ""),
        "blastx_qcovhsp": row.get("blastx_qcovhsp", ""),
        "combined_score": row.get("combined_score", ""),
        "host_source": row.get("host_source", ""),
        "plant_virus_status": row.get("plant_virus_status", "unknown"),
        "retention_reason": reason,
    }


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_ids(path, rows):
    with Path(path).open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(row["representative_query_id"] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Select generic viral candidates from virus_candidate_summary.tsv."
    )
    parser.add_argument("--summary", required=True, help="Input virus_candidate_summary.tsv")
    parser.add_argument("--out-tsv", required=True, help="Output retained candidate TSV")
    parser.add_argument("--out-ids", required=True, help="Output retained candidate ID list")
    parser.add_argument("--secondary-min-qlen", type=int, default=5000)
    parser.add_argument("--secondary-min-qcovhsp", type=float, default=70.0)
    parser.add_argument("--max-secondary", type=int, default=1)
    parser.add_argument("--secondary-backup-per-role", type=int, default=1)
    parser.add_argument("--vmr-xlsx", help="Optional ICTV VMR xlsx file for genome type annotation")
    args = parser.parse_args()

    rows = load_rows(args.summary)

    retained = []
    used_ids = set()

    primary_groups = group_rows_by_virus(rows, is_generic_primary_candidate)
    ordered_primary_groups = []
    for virus_key, virus_rows in primary_groups.items():
        best = choose_best(virus_rows, lambda row: True, rank_row_supported)
        if best is not None:
            ordered_primary_groups.append((virus_key, best, virus_rows))
    ordered_primary_groups.sort(key=lambda item: rank_row_supported(item[1]), reverse=True)

    for virus_key, best_row, virus_rows in ordered_primary_groups:
        retained.append(
            build_output_row(
                "primary",
                primary_candidate_name(best_row),
                infer_putative_virus_name(best_row),
                infer_segment_or_protein(best_row),
                best_row,
                "Top generic viral candidate for this putative virus, selected from the strongest supported rows in the summary table",
            )
        )
        used_ids.add(best_row["query_id"])

        if args.secondary_backup_per_role > 0:
            backups = choose_top_n(
                virus_rows,
                lambda row: True,
                rank_row_supported,
                args.secondary_backup_per_role,
                exclude_ids=used_ids,
            )
            for index, backup in enumerate(backups, 1):
                retained.append(
                    build_output_row(
                        "secondary",
                        backup_candidate_name(backup, index),
                        infer_putative_virus_name(backup),
                        infer_segment_or_protein(backup),
                        backup,
                        "Backup contig for the same putative virus, retained because transcriptome assembly may miss full-length viral genomes or segments",
                    )
                )
                used_ids.add(backup["query_id"])

    extra_secondary_pool = [
        row
        for row in rows
        if row.get("query_id") not in used_ids
        and virus_group_key(row) not in primary_groups
        and is_generic_secondary_candidate(row, args.secondary_min_qlen, args.secondary_min_qcovhsp)
    ]
    extra_secondary_pool.sort(key=rank_row_supported, reverse=True)

    seen_extra_viruses = set()
    extra_index = 1
    for row in extra_secondary_pool:
        key = virus_group_key(row)
        if key in seen_extra_viruses:
            continue
        retained.append(
            build_output_row(
                "secondary",
                extra_secondary_candidate_name(row, extra_index),
                infer_putative_virus_name(row),
                infer_segment_or_protein(row),
                row,
                "Additional generic viral candidate retained outside the primary set because it still shows credible viral support",
            )
        )
        seen_extra_viruses.add(key)
        extra_index += 1
        if len(seen_extra_viruses) >= args.max_secondary:
            break

    if not retained:
        generic_rows = [row for row in rows if is_generic_fallback_candidate(row)]
        generic_rows.sort(key=rank_row_supported, reverse=True)
        for index, row in enumerate(generic_rows[:10], 1):
            retained.append(
                build_output_row(
                    "secondary",
                    generic_candidate_name(row, index),
                    infer_putative_virus_name(row),
                    infer_segment_or_protein(row),
                    row,
                    "No primary candidates passed the generic rules, so this fallback viral candidate was retained from the summary table",
                )
            )

    fieldnames = [
        "priority",
        "candidate_name",
        "representative_query_id",
        "group_id",
        "putative_virus",
        "putative_segment_or_protein",
        "status",
        "best_nt_hit",
        "blastn_pident",
        "blastn_qcovs",
        "best_protein_hit",
        "blastx_pident",
        "blastx_qcovhsp",
        "combined_score",
        "genome_composition",
        "virus_genome_type_cn",
        "baltimore_group",
        "host_source",
        "plant_virus_status",
        "retention_reason",
    ]

    if args.vmr_xlsx:
        annotate_with_vmr(retained, args.vmr_xlsx)
    else:
        for row in retained:
            row["genome_composition"] = ""
            row["virus_genome_type_cn"] = ""
            row["baltimore_group"] = ""
            row["host_source"] = ""
            row["plant_virus_status"] = "unknown"

    write_tsv(args.out_tsv, retained, fieldnames)
    write_ids(args.out_ids, retained)


if __name__ == "__main__":
    main()
