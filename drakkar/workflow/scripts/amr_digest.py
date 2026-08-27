#!/usr/bin/env python3
"""Normalize AMR callers and attach geNomad mobility context.

The evidence table is deliberately lossless: AMRFinderPlus and CARD/RGI rows
remain independently identifiable.  Loci are a coordinate-based projection
for analysis, not a claim that the two databases share an ontology.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path


HIT_COLUMNS = [
    "assembly_id", "hit_id", "locus_id", "source", "source_hit_id", "gene_id",
    "contig", "start", "end", "strand", "gene_symbol", "gene_name",
    "ontology_id", "drug_class", "drug_subclass", "resistance_mechanism",
    "gene_family", "model_type", "detection_grade", "method", "mutation",
    "identity", "reference_coverage", "bitscore", "threshold", "is_partial",
    "raw_details",
]

LOCUS_COLUMNS = [
    "assembly_id", "locus_id", "contig", "start", "end", "strand",
    "primary_gene", "gene_symbols", "gene_families", "ontology_ids",
    "drug_classes", "drug_subclasses", "resistance_mechanisms", "sources",
    "source_count", "hit_count", "support_status", "concordance", "details",
]

REGION_COLUMNS = [
    "assembly_id", "region_id", "contig", "start", "end", "context_type",
    "seq_name", "length", "topology", "score", "fdr", "marker_enrichment",
    "hallmark_count", "gene_count", "conjugation_genes", "amr_accessions",
    "taxonomy", "raw_details",
]

MOBILITY_COLUMNS = [
    "assembly_id", "locus_id", "region_id", "context_type", "contig",
    "overlap_bp", "locus_overlap_fraction", "region_score",
]

DIGEST_COLUMNS = [
    *LOCUS_COLUMNS,
    "mobility_status", "mobility_region_ids", "mobility_scores",
]

DRUG_COLUMNS = [
    "assembly_id", "locus_id", "hit_id", "source", "drug_class",
    "drug_subclass", "resistance_mechanism", "gene_family", "ontology_id",
]


def clean(value):
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() in {"", "na", "nan", "none", "null"} else value


def number(value, integer=False):
    value = clean(value)
    if not value:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if not math.isfinite(parsed):
        return None
    return int(parsed) if integer else parsed


def format_number(value):
    if value is None:
        return ""
    if isinstance(value, int) or float(value).is_integer():
        return str(int(value))
    return f"{float(value):.10g}"


def canonical_key(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def row_value(row, *aliases):
    indexed = {canonical_key(key): value for key, value in row.items()}
    for alias in aliases:
        value = clean(indexed.get(canonical_key(alias)))
        if value:
            return value
    return ""


def compact_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def stable_id(prefix, *parts):
    digest = hashlib.sha256("\x1f".join(map(str, parts)).encode("utf-8")).hexdigest()[:16]
    return f"{prefix}:{digest}"


def normalize_strand(value):
    value = clean(value)
    if value in {"+", "1", "+1"}:
        return "+"
    if value in {"-", "-1"}:
        return "-"
    return ""


def normalize_interval(start, end):
    start = number(start, integer=True)
    end = number(end, integer=True)
    if start is None or end is None:
        return "", ""
    return (start, end) if start <= end else (end, start)


def read_tsv_table(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return [], []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def read_tsv(path):
    return read_tsv_table(path)[0]


def require_columns(path, fieldnames, label, *alias_groups):
    present = {canonical_key(field) for field in fieldnames}
    missing = [
        "/".join(group)
        for group in alias_groups
        if not any(canonical_key(alias) in present for alias in group)
    ]
    if missing:
        raise ValueError(
            f"{label} output {path} is missing required column(s): {', '.join(missing)}"
        )


def write_tsv(path, columns, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_amrfinder(path, assembly_id):
    hits = []
    rows, fields = read_tsv_table(path)
    require_columns(
        path, fields, "AMRFinderPlus",
        ("Contig id", "Contig"), ("Start",), ("Stop", "End"),
        ("Element symbol", "Gene symbol"), ("Type", "Element type"),
    )
    for index, native in enumerate(rows, start=1):
        element_type = row_value(native, "Type", "Element type")
        if element_type and element_type.upper() != "AMR":
            continue
        start, end = normalize_interval(
            row_value(native, "Start"), row_value(native, "Stop", "End")
        )
        method = row_value(native, "Method")
        subtype = row_value(native, "Subtype", "Element subtype")
        symbol = row_value(native, "Element symbol", "Gene symbol")
        source_hit_id = f"amrfinderplus:{index}"
        hits.append({
            "assembly_id": assembly_id,
            "hit_id": stable_id(assembly_id, "amrfinderplus", index, native),
            "source": "amrfinderplus",
            "source_hit_id": source_hit_id,
            "gene_id": row_value(native, "Protein identifier", "Protein id"),
            "contig": row_value(native, "Contig id", "Contig"),
            "start": start,
            "end": end,
            "strand": normalize_strand(row_value(native, "Strand")),
            "gene_symbol": symbol,
            "gene_name": row_value(native, "Element name", "Sequence name"),
            "ontology_id": row_value(
                native, "Hierarchy node", "HMM accession", "HMM id",
                "Closest reference accession", "Accession of closest sequence"
            ),
            "drug_class": row_value(native, "Class"),
            "drug_subclass": row_value(native, "Subclass"),
            "resistance_mechanism": "",
            "gene_family": row_value(native, "HMM description") or symbol,
            "model_type": subtype,
            "detection_grade": method,
            "method": "amrfinderplus",
            "mutation": symbol if subtype.upper().startswith("POINT") else "",
            "identity": row_value(
                native, "% Identity to reference", "% Identity to reference sequence"
            ),
            "reference_coverage": row_value(
                native, "% Coverage of reference", "% Coverage of reference sequence"
            ),
            "bitscore": "",
            "threshold": "",
            "is_partial": str("PARTIAL" in method.upper() or "INTERNAL_STOP" in method.upper()).lower(),
            "raw_details": compact_json(native),
        })
    return hits


def parse_rgi(path, assembly_id):
    hits = []
    rows, fields = read_tsv_table(path)
    require_columns(
        path, fields, "RGI",
        ("ORF_ID", "ORF ID"), ("Contig",), ("Start",), ("Stop", "End"),
        ("Best_Hit_ARO", "Best Hit ARO"), ("ARO",), ("Cut_Off", "Cut Off"),
    )
    for index, native in enumerate(rows, start=1):
        start, end = normalize_interval(
            row_value(native, "Start"), row_value(native, "Stop", "End")
        )
        cutoff = row_value(native, "Cut_Off", "Cut Off")
        source_hit_id = f"rgi:{index}"
        hits.append({
            "assembly_id": assembly_id,
            "hit_id": stable_id(assembly_id, "rgi", index, native),
            "source": "rgi",
            "source_hit_id": source_hit_id,
            "gene_id": row_value(native, "ORF_ID", "ORF ID"),
            "contig": row_value(native, "Contig"),
            "start": start,
            "end": end,
            "strand": normalize_strand(row_value(native, "Orientation", "Strand")),
            "gene_symbol": row_value(native, "Best_Hit_ARO", "Best Hit ARO"),
            "gene_name": row_value(native, "Best_Hit_ARO", "Best Hit ARO"),
            "ontology_id": row_value(native, "ARO"),
            "drug_class": row_value(native, "Drug Class"),
            "drug_subclass": row_value(native, "Antibiotic"),
            "resistance_mechanism": row_value(native, "Resistance Mechanism"),
            "gene_family": row_value(native, "AMR Gene Family"),
            "model_type": row_value(native, "Model_type", "Model type"),
            "detection_grade": cutoff,
            "method": "rgi_main",
            "mutation": row_value(native, "SNPs_in_Best_Hit_ARO", "SNPs in Best Hit ARO"),
            "identity": row_value(native, "Best_Identities", "Best Identities"),
            "reference_coverage": row_value(native, "Percentage Length of Reference Sequence"),
            "bitscore": row_value(native, "Best_Hit_bit-score", "Best Hit bit-score"),
            "threshold": row_value(native, "Pass_bit-score", "Pass bit-score"),
            "is_partial": str("partial" in row_value(native, "Note").lower()).lower(),
            "raw_details": compact_json(native),
        })
    return hits


def split_values(value):
    return [item.strip() for item in re.split(r"\s*[;,]\s*", clean(value)) if item.strip()]


def genomad_parent(gene_id):
    return re.sub(r"_\d+$", "", clean(gene_id))


def genomad_amr_by_sequence(*gene_paths):
    values = defaultdict(set)
    for path in gene_paths:
        for row in read_tsv(path):
            parent = genomad_parent(row_value(row, "gene", "seq_name"))
            for accession in split_values(row_value(row, "annotation_amr")):
                values[parent].add(accession)
    return values


def parse_coordinates(value):
    match = re.fullmatch(r"(\d+)-(\d+)", clean(value))
    if not match:
        return "", ""
    return normalize_interval(match.group(1), match.group(2))


def parse_genomad(plasmid_summary, plasmid_genes, virus_summary, virus_genes, assembly_id):
    _, plasmid_fields = read_tsv_table(plasmid_summary)
    _, plasmid_gene_fields = read_tsv_table(plasmid_genes)
    _, virus_fields = read_tsv_table(virus_summary)
    _, virus_gene_fields = read_tsv_table(virus_genes)
    require_columns(
        plasmid_summary, plasmid_fields, "geNomad plasmid summary",
        ("seq_name",), ("length",), ("plasmid_score",),
    )
    require_columns(
        virus_summary, virus_fields, "geNomad virus summary",
        ("seq_name",), ("length",), ("topology",), ("coordinates",), ("virus_score",),
    )
    require_columns(
        plasmid_genes, plasmid_gene_fields, "geNomad plasmid genes",
        ("gene",), ("annotation_amr",),
    )
    require_columns(
        virus_genes, virus_gene_fields, "geNomad virus genes",
        ("gene",), ("annotation_amr",),
    )
    amr_by_sequence = genomad_amr_by_sequence(plasmid_genes, virus_genes)
    regions = []

    for native in read_tsv(plasmid_summary):
        seq_name = row_value(native, "seq_name")
        length = number(row_value(native, "length"), integer=True)
        start, end = (1, length) if length is not None else ("", "")
        amr_accessions = set(amr_by_sequence.get(seq_name, set()))
        amr_accessions.update(split_values(row_value(native, "amr_genes")))
        regions.append({
            "assembly_id": assembly_id,
            "region_id": stable_id(assembly_id, "genomad", "plasmid", seq_name),
            "contig": seq_name,
            "start": start,
            "end": end,
            "context_type": "plasmid",
            "seq_name": seq_name,
            "length": length or "",
            "topology": row_value(native, "topology"),
            "score": row_value(native, "plasmid_score"),
            "fdr": row_value(native, "fdr"),
            "marker_enrichment": row_value(native, "marker_enrichment"),
            "hallmark_count": row_value(native, "n_hallmarks"),
            "gene_count": row_value(native, "n_genes"),
            "conjugation_genes": row_value(native, "conjugation_genes"),
            "amr_accessions": ";".join(sorted(amr_accessions)),
            "taxonomy": "",
            "raw_details": compact_json(native),
        })

    for native in read_tsv(virus_summary):
        seq_name = row_value(native, "seq_name")
        topology = row_value(native, "topology")
        length = number(row_value(native, "length"), integer=True)
        if topology.lower() == "provirus":
            context_type = "provirus"
            start, end = parse_coordinates(row_value(native, "coordinates"))
            contig = seq_name.split("|provirus_", 1)[0]
        else:
            context_type = "virus"
            start, end = (1, length) if length is not None else ("", "")
            contig = seq_name
        regions.append({
            "assembly_id": assembly_id,
            "region_id": stable_id(assembly_id, "genomad", context_type, seq_name),
            "contig": contig,
            "start": start,
            "end": end,
            "context_type": context_type,
            "seq_name": seq_name,
            "length": length or "",
            "topology": topology,
            "score": row_value(native, "virus_score"),
            "fdr": row_value(native, "fdr"),
            "marker_enrichment": row_value(native, "marker_enrichment"),
            "hallmark_count": row_value(native, "n_hallmarks"),
            "gene_count": row_value(native, "n_genes"),
            "conjugation_genes": "",
            "amr_accessions": ";".join(sorted(amr_by_sequence.get(seq_name, set()))),
            "taxonomy": row_value(native, "taxonomy"),
            "raw_details": compact_json(native),
        })
    return sorted(regions, key=lambda row: (row["contig"], int(row["start"] or 0), row["context_type"]))


def overlap_length(left_start, left_end, right_start, right_end):
    if any(value in (None, "") for value in (left_start, left_end, right_start, right_end)):
        return 0
    return max(0, min(int(left_end), int(right_end)) - max(int(left_start), int(right_start)) + 1)


def calls_match(left, right, minimum_overlap):
    # Reconciliation joins independent callers. Multiple overlapping models
    # emitted by one caller remain separate evidence loci unless a call from
    # the other source bridges them explicitly.
    if left["source"] == right["source"]:
        return False
    if not left["contig"] or left["contig"] != right["contig"]:
        return False
    if left["strand"] and right["strand"] and left["strand"] != right["strand"]:
        return False
    overlap = overlap_length(left["start"], left["end"], right["start"], right["end"])
    if overlap == 0:
        return False
    shorter = min(
        int(left["end"]) - int(left["start"]) + 1,
        int(right["end"]) - int(right["start"]) + 1,
    )
    return overlap / shorter >= minimum_overlap


def unique_values(rows, column):
    return sorted({clean(row.get(column)) for row in rows if clean(row.get(column))})


def normalized_names(rows, source):
    names = set()
    for row in rows:
        if row["source"] != source:
            continue
        for column in ("gene_symbol", "gene_family"):
            value = canonical_key(row.get(column, ""))
            if value:
                names.add(value)
    return names


def drug_terms(rows, source):
    return {
        canonical_key(term)
        for row in rows if row["source"] == source
        for term in split_values(row.get("drug_class"))
        if canonical_key(term)
    }


def concordance(rows):
    sources = {row["source"] for row in rows}
    if sources != {"amrfinderplus", "rgi"}:
        return "single_source"
    if normalized_names(rows, "amrfinderplus") & normalized_names(rows, "rgi"):
        return "gene_match"
    if drug_terms(rows, "amrfinderplus") & drug_terms(rows, "rgi"):
        return "drug_class_match"
    return "unresolved"


def reconcile_hits(hits, assembly_id, minimum_overlap):
    parent = list(range(len(hits)))

    def find(item):
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    by_contig = defaultdict(list)
    for index, hit in enumerate(hits):
        if hit["contig"] and hit["start"] != "" and hit["end"] != "":
            by_contig[hit["contig"]].append(index)
    for indexes in by_contig.values():
        for offset, left_index in enumerate(indexes):
            for right_index in indexes[offset + 1:]:
                if calls_match(hits[left_index], hits[right_index], minimum_overlap):
                    union(left_index, right_index)

    groups = defaultdict(list)
    for index, hit in enumerate(hits):
        groups[find(index)].append(hit)

    loci = []
    hit_to_locus = {}
    for rows in groups.values():
        contigs = unique_values(rows, "contig")
        starts = [int(row["start"]) for row in rows if row["start"] != ""]
        ends = [int(row["end"]) for row in rows if row["end"] != ""]
        strands = unique_values(rows, "strand")
        sources = unique_values(rows, "source")
        symbols = unique_values(rows, "gene_symbol")
        families = unique_values(rows, "gene_family")
        ontology_ids = unique_values(rows, "ontology_id")
        locus_id = stable_id(
            assembly_id,
            "locus",
            contigs[0] if len(contigs) == 1 else "",
            min(starts) if starts else "",
            max(ends) if ends else "",
            *sorted(row["hit_id"] for row in rows),
        )
        support_status = (
            "amrfinder_and_rgi"
            if set(sources) == {"amrfinderplus", "rgi"}
            else "amrfinder_only" if sources == ["amrfinderplus"]
            else "rgi_only" if sources == ["rgi"]
            else "other"
        )
        primary = next(
            (row["gene_symbol"] for row in rows if row["source"] == "amrfinderplus" and row["gene_symbol"]),
            symbols[0] if symbols else (families[0] if families else ""),
        )
        locus = {
            "assembly_id": assembly_id,
            "locus_id": locus_id,
            "contig": contigs[0] if len(contigs) == 1 else "",
            "start": min(starts) if starts else "",
            "end": max(ends) if ends else "",
            "strand": strands[0] if len(strands) == 1 else "",
            "primary_gene": primary,
            "gene_symbols": ";".join(symbols),
            "gene_families": ";".join(families),
            "ontology_ids": ";".join(ontology_ids),
            "drug_classes": ";".join(unique_values(rows, "drug_class")),
            "drug_subclasses": ";".join(unique_values(rows, "drug_subclass")),
            "resistance_mechanisms": ";".join(unique_values(rows, "resistance_mechanism")),
            "sources": ";".join(sources),
            "source_count": len(sources),
            "hit_count": len(rows),
            "support_status": support_status,
            "concordance": concordance(rows),
            "details": compact_json({"hit_ids": sorted(row["hit_id"] for row in rows)}),
        }
        loci.append(locus)
        for row in rows:
            hit_to_locus[row["hit_id"]] = locus_id
    loci.sort(key=lambda row: (row["contig"], int(row["start"] or 0), row["locus_id"]))
    return loci, hit_to_locus


def attach_mobility(loci, regions):
    links = []
    by_locus = defaultdict(list)
    for locus in loci:
        if not locus["contig"] or locus["start"] == "" or locus["end"] == "":
            continue
        locus_length = int(locus["end"]) - int(locus["start"]) + 1
        for region in regions:
            if region["contig"] != locus["contig"]:
                continue
            overlap = overlap_length(locus["start"], locus["end"], region["start"], region["end"])
            if overlap == 0:
                continue
            link = {
                "assembly_id": locus["assembly_id"],
                "locus_id": locus["locus_id"],
                "region_id": region["region_id"],
                "context_type": region["context_type"],
                "contig": locus["contig"],
                "overlap_bp": overlap,
                "locus_overlap_fraction": format_number(overlap / locus_length),
                "region_score": region["score"],
            }
            links.append(link)
            by_locus[locus["locus_id"]].append(link)

    digests = []
    for locus in loci:
        matched = by_locus.get(locus["locus_id"], [])
        context_types = sorted({row["context_type"] for row in matched})
        digests.append({
            **locus,
            "mobility_status": ";".join(context_types) if context_types else "unclassified",
            "mobility_region_ids": ";".join(sorted({row["region_id"] for row in matched})),
            "mobility_scores": ";".join(
                sorted({clean(row["region_score"]) for row in matched if clean(row["region_score"])})
            ),
        })
    return links, digests


def drug_rows(hits, hit_to_locus):
    rows = []
    seen = set()
    for hit in hits:
        classes = split_values(hit["drug_class"]) or [""]
        # Subclasses do not necessarily map one-to-one to CARD's classes.
        # Keep the native packed value rather than manufacturing a Cartesian
        # class/subclass relationship that neither caller asserted.
        subclass = clean(hit["drug_subclass"])
        for drug_class in classes:
            row = {
                "assembly_id": hit["assembly_id"],
                "locus_id": hit_to_locus[hit["hit_id"]],
                "hit_id": hit["hit_id"],
                "source": hit["source"],
                "drug_class": drug_class,
                "drug_subclass": subclass,
                "resistance_mechanism": hit["resistance_mechanism"],
                "gene_family": hit["gene_family"],
                "ontology_id": hit["ontology_id"],
            }
            signature = tuple(row[column] for column in DRUG_COLUMNS)
            if signature not in seen:
                rows.append(row)
                seen.add(signature)
    return rows


def digest(args):
    amrfinder_hits = parse_amrfinder(args.amrfinder, args.assembly_id)
    rgi_hits = parse_rgi(args.rgi, args.assembly_id)
    hits = sorted(
        [*amrfinder_hits, *rgi_hits],
        key=lambda row: (row["contig"], int(row["start"] or 0), row["source"], row["hit_id"]),
    )
    regions = parse_genomad(
        args.plasmid_summary,
        args.plasmid_genes,
        args.virus_summary,
        args.virus_genes,
        args.assembly_id,
    )
    loci, hit_to_locus = reconcile_hits(hits, args.assembly_id, args.minimum_overlap)
    for hit in hits:
        hit["locus_id"] = hit_to_locus[hit["hit_id"]]
    mobility, digests = attach_mobility(loci, regions)
    drugs = drug_rows(hits, hit_to_locus)

    write_tsv(args.hits_output, HIT_COLUMNS, hits)
    write_tsv(args.loci_output, LOCUS_COLUMNS, loci)
    write_tsv(args.regions_output, REGION_COLUMNS, regions)
    write_tsv(args.mobility_output, MOBILITY_COLUMNS, mobility)
    write_tsv(args.drugs_output, DRUG_COLUMNS, drugs)
    write_tsv(args.digest_output, DIGEST_COLUMNS, digests)

    qc = {
        "schema_version": "drakkar-amr-qc-v1",
        "assembly_id": args.assembly_id,
        "minimum_locus_overlap": args.minimum_overlap,
        "sources": [
            {
                "source": "amrfinderplus",
                "retained_hits": len(amrfinder_hits),
                "hits_without_coordinates": sum(not row["contig"] or row["start"] == "" for row in amrfinder_hits),
            },
            {
                "source": "rgi",
                "retained_hits": len(rgi_hits),
                "hits_without_coordinates": sum(not row["contig"] or row["start"] == "" for row in rgi_hits),
            },
            {"source": "genomad", "retained_regions": len(regions)},
        ],
        "loci": len(loci),
        "multi_tool_loci": sum(row["support_status"] == "amrfinder_and_rgi" for row in loci),
        "mobility_links": len(mobility),
        "mobile_loci": len({row["locus_id"] for row in mobility}),
    }
    Path(args.qc_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.qc_output).write_text(compact_json(qc) + "\n", encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Normalize AMRFinderPlus/RGI calls and attach geNomad context."
    )
    parser.add_argument("--assembly-id", required=True)
    parser.add_argument("--amrfinder", required=True)
    parser.add_argument("--rgi", required=True)
    parser.add_argument("--plasmid-summary", required=True)
    parser.add_argument("--plasmid-genes", required=True)
    parser.add_argument("--virus-summary", required=True)
    parser.add_argument("--virus-genes", required=True)
    parser.add_argument("--minimum-overlap", type=float, default=0.8)
    parser.add_argument("--hits-output", required=True)
    parser.add_argument("--loci-output", required=True)
    parser.add_argument("--regions-output", required=True)
    parser.add_argument("--mobility-output", required=True)
    parser.add_argument("--drugs-output", required=True)
    parser.add_argument("--digest-output", required=True)
    parser.add_argument("--qc-output", required=True)
    args = parser.parse_args()
    if args.minimum_overlap <= 0 or args.minimum_overlap > 1:
        parser.error("--minimum-overlap must be greater than 0 and at most 1")
    return args


if __name__ == "__main__":
    digest(parse_args())
