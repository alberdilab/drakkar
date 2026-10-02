import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml
from Bio import SearchIO

# The workflow scripts directory is not an installed package. Running this file
# directly puts it on sys.path, but importlib-based loaders (the test suite) do
# not, so add it explicitly before importing a sibling module.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from amr_columns import (
    AMRFINDER_COLUMNS,
    AMRFINDER_REQUIRED,
    RGI_COLUMNS,
    RGI_REQUIRED,
    cutoff_rank,
    field,
    indexed_row,
    method_rank,
    missing_columns,
)


DEFAULT_EVALUE_THRESHOLD = 1e-10
DEFAULT_IDENTITY_THRESHOLD = 50.0
DEFAULT_QUERY_COVERAGE_THRESHOLD = 0.5
DEFAULT_TARGET_COVERAGE_THRESHOLD = 0.5
VFDB_MAPPING_SCHEMA = "drakkar-vfdb-v2"
GENE_QC_SCHEMA = "drakkar-gene-annotation-qc-v1"

HIT_COLUMNS = [
    "gene", "source", "method", "evidence", "hit_rank", "is_primary",
    "rank_score", "rank_score_type",
    "annotation_id", "annotation", "annotation_type", "evalue", "bitscore",
    "score", "score_type", "threshold", "identity", "coverage",
    "query_coverage", "target_coverage", "confidence",
    "alignment_length", "query_start", "query_end", "target_start", "target_end",
    "model_start", "model_end", "details",
]

OUTPUT_COLUMNS = [
    "mag", "gene", "contig", "start", "end", "strand", *HIT_COLUMNS[1:],
]

NUMERIC_HIT_COLUMNS = [
    "evalue", "bitscore", "score", "threshold", "identity", "coverage",
    "query_coverage", "target_coverage", "confidence", "rank_score",
    "alignment_length", "query_start", "query_end", "target_start", "target_end",
    "model_start", "model_end",
]

SOURCE_ORDER = {
    "prodigal": 0,
    "kegg": 10,
    "pfam": 20,
    "ncbifam": 25,
    "cazy": 30,
    "vfdb": 40,
    "ncbi_amrfinder": 50,
    "card": 55,
    "signalp": 60,
    "defensefinder": 70,
    "uniprot_swissprot": 80,
}

SOURCE_RANKING = {
    "kegg": [("rank_score", False), ("score", False), ("evalue", True)],
    "pfam": [("bitscore", False), ("evalue", True)],
    "ncbifam": [
        ("rank_score", False), ("bitscore", False), ("evalue", True),
    ],
    "cazy": [("coverage", False), ("evalue", True)],
    "vfdb": [
        ("coverage", False), ("identity", False), ("bitscore", False),
        ("evalue", True),
    ],
    # AMRFinderPlus and RGI both publish an evidence tier of their own, so a
    # gene's competing calls are ranked the way each tool would rank them
    # rather than by raw alignment score.
    "ncbi_amrfinder": [
        ("rank_score", False), ("identity", False), ("coverage", False),
    ],
    "card": [
        ("rank_score", False), ("bitscore", False), ("identity", False),
    ],
    "signalp": [("confidence", False)],
    "defensefinder": [("score", False), ("evalue", True)],
    "uniprot_swissprot": [
        ("bitscore", False), ("identity", False), ("evalue", True),
    ],
    "prodigal": [("score", False)],
}


def has_content(path):
    if not path:
        return False
    path_obj = Path(path)
    return path_obj.is_file() and path_obj.stat().st_size > 0


def empty_hits():
    return pd.DataFrame(columns=HIT_COLUMNS)


def json_safe(value):
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            cleaned = json_safe(item)
            if cleaned is not None:
                result[str(key)] = cleaned
        return result
    if isinstance(value, (list, tuple, set)):
        result = []
        for item in value:
            cleaned = json_safe(item)
            if cleaned is not None:
                result.append(cleaned)
        return result
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def details_json(values):
    cleaned = json_safe(values)
    return json.dumps(cleaned or {}, sort_keys=True, separators=(",", ":"))


def first_nonempty(mapping, keys, default=""):
    for key in keys:
        value = mapping.get(key)
        cleaned = json_safe(value)
        if cleaned not in (None, ""):
            return cleaned
    return default


def negative_log10(value):
    value = pd.to_numeric(value, errors="coerce")
    if pd.isna(value):
        return pd.NA
    if value <= 0:
        return 300.0
    return -math.log10(value)


def attach_qc(
    frame,
    source,
    reported_hits,
    rejected_hits=None,
    unmapped_hits=0,
    filter_stage="drakkar",
):
    frame.attrs["annotation_qc"] = {
        "source": source,
        "reported_records": int(reported_hits),
        "retained_records": int(len(frame)),
        "rejected_records": None if rejected_hits is None else int(rejected_hits),
        "unmapped_records": int(unmapped_hits),
        "unique_entities": int(frame["gene"].nunique()) if "gene" in frame.columns else 0,
        "filter_stage": filter_stage,
    }
    return frame


def finalize_hits(
    frame,
    source,
    method,
    evidence,
    *,
    reported_hits=None,
    rejected_hits=None,
    unmapped_hits=0,
    filter_stage="drakkar",
):
    """Normalize one source and rank every retained hit without discarding it."""
    if frame is None or frame.empty:
        return attach_qc(
            empty_hits(), source, reported_hits or 0, rejected_hits,
            unmapped_hits, filter_stage,
        )

    frame = frame.copy()
    frame["source"] = source
    frame["method"] = method
    frame["evidence"] = evidence
    for column in HIT_COLUMNS:
        if column not in frame.columns:
            frame[column] = pd.NA

    frame["gene"] = frame["gene"].astype("string").str.strip()
    frame = frame[frame["gene"].notna() & (frame["gene"] != "")].copy()
    if frame.empty:
        return attach_qc(
            empty_hits(), source, reported_hits or 0, rejected_hits,
            unmapped_hits, filter_stage,
        )

    for column in NUMERIC_HIT_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")

    frame["details"] = frame["details"].fillna("{}").astype(str)
    if frame["rank_score"].isna().all():
        frame["rank_score"] = frame["score"].fillna(frame["bitscore"])
    frame["rank_score_type"] = frame["rank_score_type"].fillna(frame["score_type"])

    sort_columns = ["gene"]
    ascending = [True]
    for index, (column, lower_is_better) in enumerate(SOURCE_RANKING.get(source, [])):
        helper = f"_rank_{index}"
        missing = float("inf") if lower_is_better else float("-inf")
        frame[helper] = frame[column].fillna(missing)
        sort_columns.append(helper)
        ascending.append(lower_is_better)
    frame["_rank_annotation_id"] = frame["annotation_id"].fillna("").astype(str)
    sort_columns.append("_rank_annotation_id")
    ascending.append(True)
    frame = frame.sort_values(sort_columns, ascending=ascending, kind="stable")
    frame["hit_rank"] = frame.groupby("gene", sort=False).cumcount() + 1
    frame["is_primary"] = frame["hit_rank"] == 1
    result = frame[HIT_COLUMNS].reset_index(drop=True)
    return attach_qc(
        result,
        source,
        len(frame) if reported_hits is None else reported_hits,
        rejected_hits,
        unmapped_hits,
        filter_stage,
    )


def parse_hmmer3_tab(path):
    fields = ["accession", "bitscore", "evalue", "id", "description", "overlap_num", "region_num"]
    columns = ["gene", *fields, "bitscore_domain"]
    if not has_content(path):
        return pd.DataFrame(columns=columns)

    hits = defaultdict(list)
    query_ids = []
    with open(path) as handle:
        for queryresult in SearchIO.parse(handle, "hmmer3-tab"):
            for hit in queryresult.hits:
                query_ids.append(queryresult.id)
                for field in fields:
                    hits[field].append(getattr(hit, field, None))
                hits["bitscore_domain"].append(hit.hsps[0].bitscore if hit.hsps else None)

    if not query_ids:
        return pd.DataFrame(columns=columns)
    data = pd.DataFrame.from_dict(hits)
    data["gene"] = query_ids
    return data


NCBIFAM_METADATA_REQUIRED = {
    "ncbi_accession",
    "label",
    "sequence_cutoff",
    "domain_cutoff",
    "hmm_length",
    "family_type",
    "source",
}


def parse_hmmer3_domtbl(path):
    """Read hmmscan ``--domtblout`` without normalizing model accessions."""
    columns = [
        "model_name", "model_accession", "model_length", "gene",
        "query_accession", "query_length", "full_evalue", "full_bitscore",
        "full_bias", "domain_number", "domain_count", "conditional_evalue",
        "independent_evalue", "domain_bitscore", "domain_bias", "model_start",
        "model_end", "query_start", "query_end", "envelope_start",
        "envelope_end", "accuracy", "model_description",
    ]
    if not has_content(path):
        return pd.DataFrame(columns=columns)

    records = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split(maxsplit=22)
            if len(fields) < 22:
                raise ValueError(
                    f"Malformed NCBIfam HMMER domtblout record at line {line_number}: "
                    f"expected at least 22 fields, found {len(fields)}"
                )
            if len(fields) == 22:
                fields.append("")
            records.append(dict(zip(columns, fields)))
    frame = pd.DataFrame(records, columns=columns)
    for column in [
        "model_length", "query_length", "full_evalue", "full_bitscore",
        "full_bias", "domain_number", "domain_count", "conditional_evalue",
        "independent_evalue", "domain_bitscore", "domain_bias", "model_start",
        "model_end", "query_start", "query_end", "envelope_start",
        "envelope_end", "accuracy",
    ]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def load_ncbifam_metadata(metadata_file):
    if not has_content(metadata_file):
        raise ValueError(
            "NCBIfam hits were found but hmm_PGAP.tsv is missing or empty. "
            "Reinstall the configured NCBIfam release; Drakkar will not use an "
            "E-value fallback."
        )
    metadata = pd.read_csv(
        metadata_file, sep="\t", header=0, dtype=str, keep_default_na=False
    )
    metadata = metadata.rename(columns={column: column.lstrip("#") for column in metadata.columns})
    missing = NCBIFAM_METADATA_REQUIRED.difference(metadata.columns)
    if missing:
        raise ValueError(
            "NCBIfam metadata is missing required column(s): "
            + ", ".join(sorted(missing))
        )
    if metadata["ncbi_accession"].duplicated().any():
        duplicate = metadata.loc[
            metadata["ncbi_accession"].duplicated(keep=False), "ncbi_accession"
        ].iloc[0]
        raise ValueError(f"NCBIfam metadata contains duplicate accession: {duplicate}")
    return metadata


def ncbifam_database_provenance(metadata_file):
    """Read the managed install manifest copied into annotation QC/details."""
    metadata_path = Path(metadata_file)
    manifest_path = metadata_path.parent / "database_versions.yaml"
    installed = {}
    if manifest_path.is_file():
        try:
            installed = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
        except (OSError, yaml.YAMLError) as error:
            raise ValueError(
                f"Unable to read NCBIfam installation manifest {manifest_path}: {error}"
            ) from error
    release = str(installed.get("requested_version") or metadata_path.parent.name)
    files = installed.get("files") if isinstance(installed.get("files"), list) else []
    checksums = {
        str(record.get("filename") or Path(str(record.get("path") or "")).name): record.get("sha256")
        for record in files
        if isinstance(record, dict) and record.get("sha256")
    }
    return {
        "release": release,
        "source_version": installed.get("source_version") or f"NCBIfam/PGAP HMM release {release}",
        "sources": installed.get("sources", []),
        "checksums": checksums,
        "installation_manifest": str(manifest_path) if manifest_path.is_file() else None,
    }


def covered_fraction(intervals, length):
    length = pd.to_numeric(length, errors="coerce")
    if pd.isna(length) or length <= 0:
        return pd.NA
    normalized = sorted(
        (min(int(start), int(end)), max(int(start), int(end)))
        for start, end in intervals
        if pd.notna(start) and pd.notna(end)
    )
    if not normalized:
        return pd.NA
    covered = 0
    current_start, current_end = normalized[0]
    for start, end in normalized[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            covered += current_end - current_start + 1
            current_start, current_end = start, end
    covered += current_end - current_start + 1
    return min(1.0, covered / float(length))


def parse_ncbifam(ncbifam_file, metadata_file):
    """Normalize exact, versioned NCBIfam hits accepted by native TC1/TC2."""
    domains = parse_hmmer3_domtbl(ncbifam_file)
    metadata = load_ncbifam_metadata(metadata_file)
    provenance = ncbifam_database_provenance(metadata_file)
    if domains.empty:
        result = attach_qc(
            empty_hits(), "ncbifam", 0, 0, filter_stage="upstream_native_validated"
        )
        result.attrs["annotation_qc"].update({
            "database_release": provenance["release"],
            "database_source_version": provenance["source_version"],
            "database_checksums": details_json(provenance["checksums"]),
        })
        return result

    by_accession = {
        str(record["ncbi_accession"]): record for record in metadata.to_dict("records")
    }
    by_label = {
        str(record["label"]): record for record in metadata.to_dict("records")
        if str(record.get("label") or "")
    }
    # HMMER uses '-' for a missing accession. Installed NCBIfam profiles carry
    # ACC, but resolving via NAME keeps malformed/custom files fail-closed
    # without ever manufacturing an unversioned identifier.
    def resolve_metadata(native):
        accession = str(native.get("model_accession") or "").strip()
        if accession in {"", "-"}:
            record = by_label.get(str(native.get("model_name") or ""))
        else:
            record = by_accession.get(accession)
        if record is None:
            raise ValueError(
                "NCBIfam HMMER output references a profile absent from "
                f"hmm_PGAP.tsv: accession={accession!r}, name={native.get('model_name')!r}"
            )
        return record

    domains = domains.copy()
    domains["_metadata"] = [resolve_metadata(row) for row in domains.to_dict("records")]
    domains["_exact_accession"] = [row["ncbi_accession"] for row in domains["_metadata"]]

    rows = []
    reported_hits = 0
    for (gene, accession), group in domains.groupby(
        ["gene", "_exact_accession"], sort=False, dropna=False
    ):
        reported_hits += 1
        metadata_row = group.iloc[0]["_metadata"]
        sequence_cutoff = pd.to_numeric(metadata_row.get("sequence_cutoff"), errors="coerce")
        domain_cutoff = pd.to_numeric(metadata_row.get("domain_cutoff"), errors="coerce")
        if pd.isna(sequence_cutoff) or pd.isna(domain_cutoff):
            raise ValueError(
                f"NCBIfam profile {accession} has no complete trusted cutoff; "
                "Drakkar does not use an E-value fallback."
            )

        full_bitscore = pd.to_numeric(group.iloc[0]["full_bitscore"], errors="coerce")
        accepted_domains = group[group["domain_bitscore"] >= domain_cutoff].copy()
        if pd.isna(full_bitscore) or full_bitscore < sequence_cutoff or accepted_domains.empty:
            # hmmscan --cut_tc should already have excluded these. Rechecking
            # protects imported/raw tables and makes the contract testable.
            continue

        accepted_domains = accepted_domains.sort_values(
            ["domain_bitscore", "independent_evalue", "domain_number"],
            ascending=[False, True, True],
            kind="stable",
        )
        best = accepted_domains.iloc[0]
        query_coverage = covered_fraction(
            zip(accepted_domains["query_start"], accepted_domains["query_end"]),
            best["query_length"],
        )
        hmm_coverage = covered_fraction(
            zip(accepted_domains["model_start"], accepted_domains["model_end"]),
            best["model_length"],
        )
        coverage_values = [value for value in (query_coverage, hmm_coverage) if pd.notna(value)]

        native_domains = []
        for native in accepted_domains.to_dict("records"):
            native_domains.append({
                "accuracy": native.get("accuracy"),
                "conditional_evalue": native.get("conditional_evalue"),
                "domain_bias": native.get("domain_bias"),
                "domain_bitscore": native.get("domain_bitscore"),
                "domain_count": native.get("domain_count"),
                "domain_number": native.get("domain_number"),
                "envelope_end": native.get("envelope_end"),
                "envelope_start": native.get("envelope_start"),
                "independent_evalue": native.get("independent_evalue"),
                "model_end": native.get("model_end"),
                "model_start": native.get("model_start"),
                "query_end": native.get("query_end"),
                "query_start": native.get("query_start"),
            })

        rows.append({
            "gene": gene,
            "annotation_id": accession,
            "annotation": first_nonempty(
                metadata_row, ["product_name", "hmm_name", "name_orig", "label"]
            ),
            "annotation_type": "protein_family",
            "evalue": best["full_evalue"],
            "bitscore": full_bitscore,
            "score": full_bitscore,
            "score_type": "full_bitscore",
            "threshold": sequence_cutoff,
            "rank_score": full_bitscore - sequence_cutoff,
            "rank_score_type": "bitscore_above_ncbifam_tc1",
            "coverage": min(coverage_values) if coverage_values else pd.NA,
            "query_coverage": query_coverage,
            "target_coverage": hmm_coverage,
            "alignment_length": (
                int(best["query_end"] - best["query_start"] + 1)
                if pd.notna(best["query_start"]) and pd.notna(best["query_end"])
                else pd.NA
            ),
            "query_start": best["query_start"],
            "query_end": best["query_end"],
            "target_start": best["model_start"],
            "target_end": best["model_end"],
            "model_start": best["model_start"],
            "model_end": best["model_end"],
            "details": details_json({
                "acceptance_rule": "hmmer_cut_tc_sequence_and_domain",
                "acceptance_status": "accepted",
                "domain_trusted_cutoff": domain_cutoff,
                "database_checksums": provenance["checksums"],
                "family_type": metadata_row.get("family_type"),
                "for_naming": metadata_row.get("for_naming"),
                "for_structural_annotation": metadata_row.get("for_structural_annotation"),
                "model_description": first_nonempty(
                    metadata_row, ["product_name", "hmm_name", "name_orig"]
                ),
                "native_domains": native_domains,
                "native_hmm": {
                    "accession": accession,
                    "description": best.get("model_description"),
                    "length": best.get("model_length"),
                    "name": best.get("model_name"),
                },
                "profile_grade": first_nonempty(metadata_row, ["profile_grade", "grade"]),
                "profile_metadata": metadata_row,
                "profile_source": metadata_row.get("source"),
                "profile_type": metadata_row.get("family_type"),
                "sequence_trusted_cutoff": sequence_cutoff,
                "source_release": provenance["release"],
                "source_urls": provenance["sources"],
                "source_version": provenance["source_version"],
                "threshold_type": "trusted_cutoff",
            }),
        })

    result = finalize_hits(
        pd.DataFrame(rows),
        "ncbifam",
        "hmmer",
        "sequence_homology",
        reported_hits=reported_hits,
        rejected_hits=reported_hits - len(rows),
        filter_stage="upstream_native_validated",
    )
    result.attrs["annotation_qc"].update({
        "database_release": provenance["release"],
        "database_source_version": provenance["source_version"],
        "database_checksums": details_json(provenance["checksums"]),
    })
    return result


def load_kofam_thresholds(kolist_file):
    columns = ["kegg", "threshold", "score_type"]
    if not has_content(kolist_file):
        return pd.DataFrame(columns=columns)
    kolist = pd.read_csv(kolist_file, sep="\t", header=0, usecols=["knum", "threshold", "score_type"])
    kolist = kolist.rename(columns={"knum": "kegg"})
    kolist["threshold"] = pd.to_numeric(kolist["threshold"], errors="coerce")
    return kolist[columns]


def load_kegg_hierarchy(keggdb_file):
    if not has_content(keggdb_file):
        return pd.DataFrame(columns=["kegg", "ec"])
    with open(keggdb_file) as handle:
        kegg_json = json.load(handle)
    rows = []
    for main in kegg_json.get("children", []):
        for broad in main.get("children", []):
            for sub in broad.get("children", []):
                for gene_node in sub.get("children", []):
                    name = gene_node.get("name", "")
                    kegg = name.split(" ")[0] if name else ""
                    description_and_ec = " ".join(name.split(" ")[1:])
                    description_parts = description_and_ec.split(" [")
                    ec = description_parts[1][:-1] if len(description_parts) > 1 else ""
                    rows.append((kegg, ec.replace("EC:", "")))

    # A KO can occur under several branches of ko00001. Joining that hierarchy
    # directly to HMMER hits would therefore turn one biological hit into many
    # output rows. Collapse the hierarchy to one record per KO while retaining
    # every distinct EC association in source order.
    hierarchy = pd.DataFrame(rows, columns=["kegg", "ec"])
    collapsed = []
    for kegg, group in hierarchy.groupby("kegg", sort=False):
        ecs = []
        seen = set()
        for value in group["ec"]:
            for ec in re.split(r"[\s;,]+", str(value).strip()):
                if ec and ec not in seen:
                    seen.add(ec)
                    ecs.append(ec)
        collapsed.append((kegg, ";".join(ecs)))
    return pd.DataFrame(collapsed, columns=["kegg", "ec"])


def parse_kegg(kegg_file, keggdb_file, kolist_file, evalue_threshold):
    hits = parse_hmmer3_tab(kegg_file)
    if hits.empty:
        return attach_qc(empty_hits(), "kegg", 0, 0)
    reported_hits = len(hits)
    for column in ["evalue", "bitscore", "bitscore_domain"]:
        hits[column] = pd.to_numeric(hits[column], errors="coerce")
    hits = hits.rename(columns={"id": "kegg"})

    thresholds = load_kofam_thresholds(kolist_file)
    if thresholds.empty:
        raise ValueError(
            "KOfam hits were found but the native ko_list cutoff table is missing or empty. "
            "Reinstall the configured KEGG/KOfam release instead of applying a global "
            "e-value fallback."
        )
    hits = pd.merge(hits, thresholds, on="kegg", how="left", validate="many_to_one")
    domain_scored = hits["score_type"] == "domain"
    selected_bitscore = hits["bitscore"].where(~domain_scored, hits["bitscore_domain"])
    has_threshold = hits["threshold"].notna()
    passes_native = has_threshold & (selected_bitscore >= hits["threshold"])
    passes_fallback = ~has_threshold & (hits["evalue"] <= evalue_threshold)
    hits["filter_score"] = selected_bitscore.where(has_threshold, hits["evalue"])
    hits["filter_score_type"] = hits["score_type"].map(
        {"full": "full_bitscore", "domain": "domain_bitscore"}
    ).where(has_threshold, "evalue")
    hits["threshold"] = hits["threshold"].where(has_threshold, evalue_threshold)
    hits["rank_score"] = (selected_bitscore - hits["threshold"]).where(
        has_threshold,
        hits["evalue"].map(negative_log10),
    )
    hits["rank_score_type"] = "bitscore_above_kofam_cutoff"
    hits.loc[~has_threshold, "rank_score_type"] = "negative_log10_evalue"
    hits = hits[passes_native | passes_fallback].copy()
    if hits.empty:
        return attach_qc(empty_hits(), "kegg", reported_hits, reported_hits)

    hierarchy = load_kegg_hierarchy(keggdb_file)
    hits = (
        pd.merge(hits, hierarchy, on="kegg", how="left", validate="many_to_one")
        if not hierarchy.empty
        else hits.assign(ec=pd.NA)
    )
    rows = []
    for _, row in hits.iterrows():
        rows.append({
            "gene": row["gene"],
            "annotation_id": row["kegg"],
            "annotation": first_nonempty(row, ["ec", "description"]),
            "annotation_type": "ko",
            "evalue": row["evalue"],
            "bitscore": row["bitscore"],
            "score": row["filter_score"],
            "score_type": row["filter_score_type"],
            "threshold": row["threshold"],
            "rank_score": row["rank_score"],
            "rank_score_type": row["rank_score_type"],
            "details": details_json({
                "accession": row.get("accession"),
                "domain_bitscore": row.get("bitscore_domain"),
                "ec": row.get("ec"),
                "hmm_description": row.get("description"),
                "overlap_num": row.get("overlap_num"),
                "region_num": row.get("region_num"),
            }),
        })
    unmapped = sum(json_safe(row.get("ec")) in (None, "") for _, row in hits.iterrows())
    return finalize_hits(
        pd.DataFrame(rows),
        "kegg",
        "hmmscan",
        "sequence_homology",
        reported_hits=reported_hits,
        rejected_hits=reported_hits - len(hits),
        unmapped_hits=unmapped,
    )


def load_pfam_ec_associations(ec_file):
    if not has_content(ec_file):
        return {}
    associations = pd.read_csv(ec_file, sep="\t", comment="#", header=0)
    associations = associations[associations["Type"] == "GOLD"].copy()
    associations["Pfam-Domain"] = associations["Pfam-Domain"].astype("string").str.split(".").str[0]
    associations["Confidence-Score"] = pd.to_numeric(associations["Confidence-Score"], errors="coerce")
    by_pfam = {}
    for pfam, group in associations.groupby("Pfam-Domain", sort=False):
        by_pfam[str(pfam)] = [
            {"ec": row.get("EC-Number"), "confidence": row.get("Confidence-Score"), "type": row.get("Type")}
            for row in group.to_dict("records")
        ]
    return by_pfam


def parse_pfam(pfam_file, ec_file):
    hits = parse_hmmer3_tab(pfam_file)
    if hits.empty:
        return attach_qc(empty_hits(), "pfam", 0, None, filter_stage="upstream_native")
    reported_hits = len(hits)
    for column in ["evalue", "bitscore", "bitscore_domain"]:
        hits[column] = pd.to_numeric(hits[column], errors="coerce")
    hits["pfam"] = hits["accession"].astype("string").str.split(".").str[0]
    associations = load_pfam_ec_associations(ec_file)
    rows = []
    for _, row in hits.iterrows():
        pfam = str(row["pfam"])
        rows.append({
            "gene": row["gene"],
            "annotation_id": pfam,
            "annotation": first_nonempty(row, ["description", "id"]),
            "annotation_type": "protein_family",
            "evalue": row["evalue"],
            "bitscore": row["bitscore"],
            "score": row["bitscore"],
            "score_type": "full_bitscore",
            "rank_score": row["bitscore"],
            "rank_score_type": "full_bitscore",
            "details": details_json({
                "domain_bitscore": row.get("bitscore_domain"),
                "ec_associations": associations.get(pfam, []),
                "model_name": row.get("id"),
                "original_accession": row.get("accession"),
                "overlap_num": row.get("overlap_num"),
                "region_num": row.get("region_num"),
                "threshold_rule": "pfam_gathering",
            }),
        })
    return finalize_hits(
        pd.DataFrame(rows),
        "pfam",
        "hmmscan",
        "sequence_homology",
        reported_hits=reported_hits,
        filter_stage="upstream_native",
    )


def parse_cazy(cazy_file):
    required = {
        "HMM Name", "HMM Length", "Target Name", "Target Length", "i-Evalue",
        "HMM From", "HMM To", "Target From", "Target To", "Coverage", "HMM File Name",
    }
    if not has_content(cazy_file):
        return attach_qc(empty_hits(), "cazy", 0, None, filter_stage="upstream_native")
    try:
        hits = pd.read_csv(cazy_file, sep="\t", comment="#", header=0)
    except pd.errors.EmptyDataError as error:
        raise ValueError(f"dbCAN CAZy output is missing required columns: {', '.join(sorted(required))}") from error
    missing = required.difference(hits.columns)
    if missing:
        raise ValueError(f"dbCAN CAZy output is missing required columns: {', '.join(sorted(missing))}")
    if hits.empty:
        return attach_qc(empty_hits(), "cazy", 0, None, filter_stage="upstream_native")
    reported_hits = len(hits)

    rows = []
    for _, row in hits.iterrows():
        family = re.sub(r"\.hmm$", "", Path(str(row["HMM Name"])).name)
        rows.append({
            "gene": row["Target Name"],
            "annotation_id": family,
            "annotation": family,
            "annotation_type": "cazy_family",
            "evalue": row["i-Evalue"],
            "score": row["Coverage"],
            "score_type": "hmm_coverage",
            "rank_score": row["Coverage"],
            "rank_score_type": "hmm_coverage",
            "threshold": 0.35,
            "coverage": row["Coverage"],
            "query_start": row["Target From"],
            "query_end": row["Target To"],
            "model_start": row["HMM From"],
            "model_end": row["HMM To"],
            "details": details_json({
                "coverage_threshold": 0.35,
                "evalue_threshold": 1e-15,
                "hmm_file": row["HMM File Name"],
                "hmm_length": row["HMM Length"],
                "target_length": row["Target Length"],
                "threshold_rule": "dbcan_native",
            }),
        })
    return finalize_hits(
        pd.DataFrame(rows),
        "cazy",
        "run_dbcan_hmm",
        "sequence_homology",
        reported_hits=reported_hits,
        filter_stage="upstream_native",
    )


def mapping_records(path, key_column, rename=None):
    if not has_content(path):
        return {}
    # AMRFinder's first column is literally named ``#hmm_accession``; treating
    # '#' as a comment marker would silently discard the real header.
    mapping = pd.read_csv(path, sep="\t", header=0)
    if rename:
        mapping = mapping.rename(columns=rename)
    if key_column not in mapping.columns:
        return {}
    grouped = {}
    for key, group in mapping.groupby(key_column, dropna=False, sort=False):
        if pd.isna(key):
            continue
        grouped[str(key)] = group.to_dict("records")
    return grouped


def load_vfdb_mapping(path):
    """Load only VFDB mappings generated with the corrected v2 parser."""
    if not has_content(path):
        raise ValueError(
            "VFDB hits were found but the VFDB mapping table is missing or empty. "
            "Install a fresh mapping with `drakkar database vfdb --directory <vfdb-root> --set-default`."
        )

    mapping = pd.read_csv(path, sep="\t", header=0)
    required = {"entry", "vf", "vfc", "vf_type", "mapping_schema"}
    missing = required.difference(mapping.columns)
    schemas = set(mapping.get("mapping_schema", pd.Series(dtype="string")).dropna().astype(str))
    if missing or schemas != {VFDB_MAPPING_SCHEMA}:
        reason = f"missing columns: {', '.join(sorted(missing))}" if missing else f"schema values: {sorted(schemas)}"
        raise ValueError(
            "VFDB mapping is incompatible with Drakkar 2.0 "
            f"({reason}). Rebuild it with `drakkar database vfdb "
            "--directory <vfdb-root> --set-default`; legacy mappings can label "
            "the organism as the virulence-factor type."
        )

    grouped = {}
    for key, group in mapping.groupby("entry", dropna=False, sort=False):
        if pd.isna(key):
            continue
        grouped[str(key)] = group.to_dict("records")
    return grouped


def read_native_tsv(path):
    """Read one native tool table, returning its records and column names."""
    if not has_content(path):
        return [], []
    frame = pd.read_csv(path, sep="\t", header=0, dtype=str, keep_default_na=False)
    return frame.to_dict("records"), list(frame.columns)


def first_token(value):
    """Take the identifier out of a header-derived field.

    Prodigal protein headers carry coordinates after the id, and RGI echoes the
    header it was given, so both tools can hand back more than the bare gene id.
    """
    return str(value).split(None, 1)[0] if str(value).strip() else ""


def parse_amr(amr_file):
    """Normalize AMRFinderPlus's native report into gene-level hits.

    Acceptance is entirely AMRFinderPlus's: its BLASTP arm applies per-gene
    curated identity and coverage cutoffs and its HMM arm applies the NCBIfam
    trusted cutoffs, so Drakkar adds no threshold of its own. Only the element
    type is filtered, keeping this source to AMR and leaving the stress and
    virulence "plus" genes to their own sources.
    """
    rows_native, fieldnames = read_native_tsv(amr_file)
    if not rows_native:
        return attach_qc(
            empty_hits(), "ncbi_amrfinder", 0, None, filter_stage="upstream_native"
        )
    missing = missing_columns(fieldnames, AMRFINDER_REQUIRED)
    if missing:
        raise ValueError(
            "AMRFinderPlus output is missing required column(s): "
            f"{', '.join(missing)}. Regenerate it with the AMRFinderPlus "
            "version pinned in workflow/envs/amr_amrfinder.yaml."
        )

    reported_hits = len(rows_native)
    rows = []
    for native in rows_native:
        indexed = indexed_row(native)
        element_type = field(indexed, AMRFINDER_COLUMNS["type"])
        if element_type and element_type.upper() != "AMR":
            continue
        method = field(indexed, AMRFINDER_COLUMNS["method"])
        identity = pd.to_numeric(
            field(indexed, AMRFINDER_COLUMNS["identity"]) or None, errors="coerce"
        )
        reference_coverage = pd.to_numeric(
            field(indexed, AMRFINDER_COLUMNS["reference_coverage"]) or None,
            errors="coerce",
        )
        # AMRFinderPlus reports coverage as a percentage; the table stores
        # coverage as a fraction so that it is comparable across sources.
        coverage = reference_coverage / 100.0 if pd.notna(reference_coverage) else pd.NA
        rows.append({
            "gene": first_token(field(indexed, AMRFINDER_COLUMNS["gene"])),
            "annotation_id": field(indexed, AMRFINDER_COLUMNS["symbol"]),
            "annotation": field(
                indexed, AMRFINDER_COLUMNS["name"],
                field(indexed, AMRFINDER_COLUMNS["symbol"]),
            ),
            "annotation_type": field(indexed, AMRFINDER_COLUMNS["drug_class"]),
            "identity": identity,
            "coverage": coverage,
            "target_coverage": coverage,
            "alignment_length": pd.to_numeric(
                field(indexed, AMRFINDER_COLUMNS["alignment_length"]) or None,
                errors="coerce",
            ),
            "score": identity,
            "score_type": "percent_identity_to_reference",
            "rank_score": method_rank(method),
            "rank_score_type": "amrfinder_method_rank",
            "details": details_json({
                "drug_class": field(indexed, AMRFINDER_COLUMNS["drug_class"]),
                "drug_subclass": field(indexed, AMRFINDER_COLUMNS["drug_subclass"]),
                "element_subtype": field(indexed, AMRFINDER_COLUMNS["subtype"]),
                "element_type": element_type,
                "hierarchy_node": field(indexed, AMRFINDER_COLUMNS["hierarchy_node"]),
                "hmm_accession": field(indexed, AMRFINDER_COLUMNS["hmm_accession"]),
                "hmm_description": field(indexed, AMRFINDER_COLUMNS["hmm_description"]),
                "method": method,
                "native": native,
                "reference_accession": field(
                    indexed, AMRFINDER_COLUMNS["reference_accession"]
                ),
                "reference_coverage_percent": field(
                    indexed, AMRFINDER_COLUMNS["reference_coverage"]
                ),
                "reference_name": field(indexed, AMRFINDER_COLUMNS["reference_name"]),
                "threshold_rule": "amrfinderplus_curated",
            }),
        })
    if not rows:
        return attach_qc(
            empty_hits(), "ncbi_amrfinder", reported_hits, reported_hits,
            filter_stage="upstream_native",
        )
    return finalize_hits(
        pd.DataFrame(rows),
        "ncbi_amrfinder",
        "amrfinderplus",
        "sequence_homology",
        reported_hits=reported_hits,
        rejected_hits=reported_hits - len(rows),
        filter_stage="upstream_native",
    )


def parse_card(card_file):
    """Normalize CARD/RGI's native report into gene-level hits.

    RGI decides acceptance with its own per-model curated bit score cutoffs and
    reports the tier it used in ``Cut_Off``. Without ``--include_loose`` it
    emits only Perfect and Strict calls, so Drakkar adds no further filter and
    keeps RGI's own cutoff in the ``threshold`` column.
    """
    rows_native, fieldnames = read_native_tsv(card_file)
    if not rows_native:
        return attach_qc(empty_hits(), "card", 0, None, filter_stage="upstream_native")
    missing = missing_columns(fieldnames, RGI_REQUIRED)
    if missing:
        raise ValueError(
            "CARD/RGI output is missing required column(s): "
            f"{', '.join(missing)}. Regenerate it with the RGI version pinned "
            "in workflow/envs/amr_rgi.yaml, run in protein mode."
        )

    reported_hits = len(rows_native)
    rows = []
    for native in rows_native:
        indexed = indexed_row(native)
        cutoff = field(indexed, RGI_COLUMNS["cutoff"])
        reference_coverage = pd.to_numeric(
            field(indexed, RGI_COLUMNS["reference_coverage"]) or None, errors="coerce"
        )
        coverage = reference_coverage / 100.0 if pd.notna(reference_coverage) else pd.NA
        bitscore = pd.to_numeric(
            field(indexed, RGI_COLUMNS["bitscore"]) or None, errors="coerce"
        )
        rows.append({
            "gene": first_token(field(indexed, RGI_COLUMNS["gene"])),
            "annotation_id": field(
                indexed, RGI_COLUMNS["aro"], field(indexed, RGI_COLUMNS["aro_name"])
            ),
            "annotation": field(indexed, RGI_COLUMNS["aro_name"]),
            "annotation_type": field(indexed, RGI_COLUMNS["drug_class"]),
            "identity": pd.to_numeric(
                field(indexed, RGI_COLUMNS["identity"]) or None, errors="coerce"
            ),
            "coverage": coverage,
            "target_coverage": coverage,
            "bitscore": bitscore,
            "score": bitscore,
            "score_type": "bitscore",
            "threshold": pd.to_numeric(
                field(indexed, RGI_COLUMNS["threshold"]) or None, errors="coerce"
            ),
            "rank_score": cutoff_rank(cutoff),
            "rank_score_type": "rgi_cutoff_rank",
            "details": details_json({
                "amr_gene_family": field(indexed, RGI_COLUMNS["gene_family"]),
                "aro": field(indexed, RGI_COLUMNS["aro"]),
                "cut_off": cutoff,
                "drug_class": field(indexed, RGI_COLUMNS["drug_class"]),
                "drug_subclass": field(indexed, RGI_COLUMNS["drug_subclass"]),
                "model_type": field(indexed, RGI_COLUMNS["model_type"]),
                "native": native,
                "note": field(indexed, RGI_COLUMNS["note"]),
                "reference_coverage_percent": field(
                    indexed, RGI_COLUMNS["reference_coverage"]
                ),
                "resistance_mechanism": field(
                    indexed, RGI_COLUMNS["resistance_mechanism"]
                ),
                "snps": field(indexed, RGI_COLUMNS["snps"]),
                "threshold_rule": "rgi_curated_bitscore",
            }),
        })
    return finalize_hits(
        pd.DataFrame(rows),
        "card",
        "rgi_main",
        "sequence_homology",
        reported_hits=reported_hits,
        filter_stage="upstream_native",
    )


def parse_vfdb(
    vf_file,
    vfdb_file,
    evalue_threshold,
    identity_threshold,
    query_coverage_threshold=DEFAULT_QUERY_COVERAGE_THRESHOLD,
    target_coverage_threshold=DEFAULT_TARGET_COVERAGE_THRESHOLD,
):
    if not has_content(vf_file):
        return attach_qc(empty_hits(), "vfdb", 0, 0)
    hits = pd.read_csv(
        vf_file,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "gene", "entry", "identity", "alignment_length", "mismatches", "gaps",
            "query_start", "query_end", "target_start", "target_end", "evalue", "bitscore",
            "query_length", "target_length", "query_coverage", "target_coverage",
        ],
    )
    if hits.empty:
        return attach_qc(empty_hits(), "vfdb", 0, 0)
    reported_hits = len(hits)
    for column in [
        "identity", "alignment_length", "query_start", "query_end", "target_start",
        "target_end", "evalue", "bitscore", "query_length", "target_length",
        "query_coverage", "target_coverage",
    ]:
        hits[column] = pd.to_numeric(hits[column], errors="coerce")
    coverage_columns = ["query_length", "target_length", "query_coverage", "target_coverage"]
    missing_coverage = [column for column in coverage_columns if hits[column].isna().all()]
    if missing_coverage:
        raise ValueError(
            "VFDB/MMseqs output does not contain the required length and coverage fields "
            f"({', '.join(missing_coverage)}). Regenerate it with Drakkar 2.0's "
            "query/target coverage format."
        )
    hits = hits[
        (hits["evalue"] <= evalue_threshold)
        & (hits["identity"] >= identity_threshold)
        & (hits["query_coverage"] >= query_coverage_threshold)
        & (hits["target_coverage"] >= target_coverage_threshold)
    ].copy()
    if hits.empty:
        return attach_qc(empty_hits(), "vfdb", reported_hits, reported_hits)

    mappings = load_vfdb_mapping(vfdb_file)
    hits["coverage"] = hits[["query_coverage", "target_coverage"]].min(axis=1)
    rows = []
    for _, row in hits.iterrows():
        mapped = mappings.get(str(row["entry"]), [])
        primary = mapped[0] if mapped else {}
        rows.append({
            "gene": row["gene"],
            "annotation_id": row["entry"],
            "annotation": primary.get("vf", ""),
            "annotation_type": primary.get("vf_type", ""),
            "evalue": row["evalue"],
            "bitscore": row["bitscore"],
            "score": row["bitscore"],
            "score_type": "bitscore",
            "rank_score": row["coverage"],
            "rank_score_type": "minimum_query_target_coverage",
            "identity": row["identity"],
            "coverage": row["coverage"],
            "query_coverage": row["query_coverage"],
            "target_coverage": row["target_coverage"],
            "alignment_length": row["alignment_length"],
            "query_start": row["query_start"],
            "query_end": row["query_end"],
            "target_start": row["target_start"],
            "target_end": row["target_end"],
            "details": details_json({
                "evalue_threshold": evalue_threshold,
                "gaps": row.get("gaps"),
                "identity_scale": "percent",
                "identity_threshold": identity_threshold,
                "query_coverage_threshold": query_coverage_threshold,
                "query_length": row.get("query_length"),
                "target_coverage_threshold": target_coverage_threshold,
                "target_length": row.get("target_length"),
                "mappings": mapped,
                "mismatches": row.get("mismatches"),
                "vfc": primary.get("vfc"),
            }),
        })
    unmapped = sum(not mappings.get(str(entry)) for entry in hits["entry"])
    return finalize_hits(
        pd.DataFrame(rows),
        "vfdb",
        "mmseqs_easy_search",
        "sequence_homology",
        reported_hits=reported_hits,
        rejected_hits=reported_hits - len(hits),
        unmapped_hits=unmapped,
    )


def parse_signalp(signalp_file):
    if not has_content(signalp_file):
        return attach_qc(empty_hits(), "signalp", 0, None, filter_stage="upstream_native")
    hits = pd.read_csv(signalp_file, sep="\t", comment="#", header=None, names=["gene", "signalp", "confidence"])
    if hits.empty:
        return attach_qc(empty_hits(), "signalp", 0, None, filter_stage="upstream_native")
    reported_hits = len(hits)
    rows = []
    for _, row in hits.iterrows():
        rows.append({
            "gene": row["gene"],
            "annotation_id": row["signalp"],
            "annotation": row["signalp"],
            "annotation_type": "signal_peptide",
            "score": row["confidence"],
            "score_type": "confidence",
            "rank_score": row["confidence"],
            "rank_score_type": "confidence",
            "confidence": row["confidence"],
            "details": "{}",
        })
    return finalize_hits(
        pd.DataFrame(rows),
        "signalp",
        "signalp6",
        "protein_feature_prediction",
        reported_hits=reported_hits,
        filter_stage="upstream_native",
    )


def parse_defensefinder(defense_file):
    if not has_content(defense_file):
        return attach_qc(
            empty_hits(), "defensefinder", 0, None, filter_stage="upstream_native"
        )
    hits = pd.read_csv(defense_file, sep="\t")
    if hits.empty or "hit_id" not in hits.columns:
        return attach_qc(
            empty_hits(), "defensefinder", 0, None, filter_stage="upstream_native"
        )
    reported_hits = len(hits)
    rows = []
    for native in hits.to_dict("records"):
        rows.append({
            "gene": native.get("hit_id"),
            "annotation_id": first_nonempty(native, ["gene_name", "profile_name", "hit_id"]),
            "annotation": first_nonempty(native, ["type", "subtype", "gene_name"]),
            "annotation_type": first_nonempty(native, ["activity"], "defense_system_gene"),
            "evalue": first_nonempty(native, ["i_evalue", "evalue"], None),
            "bitscore": first_nonempty(native, ["bitscore", "hit_score"], None),
            "score": first_nonempty(native, ["hit_score", "bitscore", "score"], None),
            "score_type": "defensefinder_score",
            "rank_score": first_nonempty(native, ["hit_score", "bitscore", "score"], None),
            "rank_score_type": "defensefinder_score",
            "details": details_json(native),
        })
    return finalize_hits(
        pd.DataFrame(rows),
        "defensefinder",
        "defense_finder",
        "defense_system_prediction",
        reported_hits=reported_hits,
        filter_stage="upstream_native",
    )


def uniprot_accession_from_target(target):
    match = re.match(r"AF-([^-]+)-F\d+", str(target))
    return match.group(1) if match else str(target)


def parse_foldseek(foldseek_file, mapdb_file, evalue_threshold):
    if not has_content(foldseek_file):
        return attach_qc(empty_hits(), "uniprot_swissprot", 0, 0)
    hits = pd.read_csv(
        foldseek_file,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "gene", "target", "identity", "alignment_length", "mismatches", "gaps",
            "query_start", "query_end", "target_start", "target_end", "evalue", "bitscore",
        ],
    )
    if hits.empty:
        return attach_qc(empty_hits(), "uniprot_swissprot", 0, 0)
    reported_hits = len(hits)
    for column in [
        "identity", "alignment_length", "query_start", "query_end", "target_start",
        "target_end", "evalue", "bitscore",
    ]:
        hits[column] = pd.to_numeric(hits[column], errors="coerce")
    hits = hits[hits["evalue"] <= evalue_threshold].copy()
    if hits.empty:
        return attach_qc(
            empty_hits(), "uniprot_swissprot", reported_hits, reported_hits
        )

    mappings = mapping_records(mapdb_file, "accession")
    rows = []
    for _, row in hits.iterrows():
        accession = uniprot_accession_from_target(row["target"])
        mapped = mappings.get(str(accession), [])
        primary = mapped[0] if mapped else {}
        mapped_label = ";".join(
            f"{key}={primary[key]}" for key in ("kegg", "ec", "pfam")
            if json_safe(primary.get(key)) not in (None, "")
        )
        rows.append({
            "gene": row["gene"],
            "annotation_id": accession,
            "annotation": mapped_label,
            "annotation_type": "structure_match",
            "evalue": row["evalue"],
            "bitscore": row["bitscore"],
            "score": row["bitscore"],
            "score_type": "bitscore",
            "rank_score": row["bitscore"],
            "rank_score_type": "bitscore",
            "identity": row["identity"],
            "alignment_length": row["alignment_length"],
            "query_start": row["query_start"],
            "query_end": row["query_end"],
            "target_start": row["target_start"],
            "target_end": row["target_end"],
            "details": details_json({
                "database": "AlphaFold/Swiss-Prot",
                "evalue_threshold": evalue_threshold,
                "gaps": row.get("gaps"),
                "identity_scale": "fraction",
                "mappings": mapped,
                "mismatches": row.get("mismatches"),
                "target": row["target"],
            }),
        })
    unmapped = sum(
        not mappings.get(str(uniprot_accession_from_target(target)))
        for target in hits["target"]
    )
    return finalize_hits(
        pd.DataFrame(rows),
        "uniprot_swissprot",
        "foldseek_prostt5",
        "structure_homology",
        reported_hits=reported_hits,
        rejected_hits=reported_hits - len(hits),
        unmapped_hits=unmapped,
    )


def gene_id_from_gff(row):
    match = re.search(r"(?:^|;)ID=([^;]+)", str(row["attributes"]))
    identifier = match.group(1) if match else str(row["attributes"]).split(";")[0].replace("ID=", "")
    return f"{row['contig']}_{identifier.split('_')[-1]}"


def parse_gene_calls(gff_file):
    columns = ["contig", "gff_source", "feature_type", "start", "end", "gff_score", "strand", "phase", "attributes"]
    try:
        genes = pd.read_csv(gff_file, sep="\t", comment="#", header=None, names=columns)
    except pd.errors.EmptyDataError:
        genes = pd.DataFrame(columns=columns)
    if genes.empty:
        calls = attach_qc(empty_hits(), "prodigal", 0, 0, filter_stage="gene_prediction")
        return pd.DataFrame(columns=["gene", "contig", "start", "end", "strand"]), calls

    genes["gene"] = genes.apply(gene_id_from_gff, axis=1)
    duplicate_ids = sorted(genes.loc[genes["gene"].duplicated(keep=False), "gene"].astype(str).unique())
    if duplicate_ids:
        preview = ", ".join(duplicate_ids[:10])
        remainder = f" (and {len(duplicate_ids) - 10} more)" if len(duplicate_ids) > 10 else ""
        raise ValueError(
            "Prodigal GFF contains duplicate derived gene IDs: "
            f"{preview}{remainder}. Gene IDs must be unique within each MAG."
        )

    metadata = genes[["gene", "contig", "start", "end", "strand"]].copy()
    rows = []
    for _, row in genes.iterrows():
        rows.append({
            "gene": row["gene"],
            "annotation_id": row["feature_type"],
            "annotation": "",
            "annotation_type": "gene_call",
            "score": row["gff_score"],
            "score_type": "gff_score",
            "rank_score": row["gff_score"],
            "rank_score_type": "gff_score",
            "details": details_json({
                "attributes": row["attributes"],
                "feature_type": row["feature_type"],
                "gff_source": row["gff_source"],
                "phase": row["phase"],
            }),
        })
    calls = finalize_hits(
        pd.DataFrame(rows),
        "prodigal",
        "prodigal",
        "gene_prediction",
        reported_hits=len(rows),
        rejected_hits=0,
        filter_stage="gene_prediction",
    )
    return metadata, calls


def validate_hit_gene_ids(hits, genes, mag):
    """Require every functional hit to resolve to exactly one Prodigal gene."""
    if hits.empty:
        return

    known_genes = set(genes["gene"].dropna().astype(str))
    functional = hits[hits["source"] != "prodigal"].copy()
    unknown = functional[~functional["gene"].astype(str).isin(known_genes)]
    if unknown.empty:
        return

    pairs = (
        unknown[["source", "gene"]]
        .drop_duplicates()
        .sort_values(["source", "gene"], kind="stable")
    )
    labels = [f"{row.source}:{row.gene}" for row in pairs.head(10).itertuples(index=False)]
    remainder = f" (and {len(pairs) - 10} more)" if len(pairs) > 10 else ""
    raise ValueError(
        f"Annotation hits for MAG {mag!r} do not match any Prodigal gene: "
        f"{', '.join(labels)}{remainder}. Check that every annotation source "
        "was generated from this MAG's Prodigal protein FASTA."
    )


GENE_SOURCE_ALIASES = {
    "vfdb": "virulence",
    "foldseek": "structure",
    "rgi": "card",
}

GENE_SOURCE_FACTORIES = {
    "kegg",
    "cazy",
    "pfam",
    "ncbifam",
    "virulence",
    "amr",
    "card",
    "signalp",
    "defense",
    "structure",
}


def normalize_enabled_sources(enabled_sources):
    if enabled_sources is None:
        # NCBIfam is intentionally opt-in even for direct script/API callers;
        # the workflow always passes an explicit source list.
        return set(GENE_SOURCE_FACTORIES).difference({"ncbifam"})
    if isinstance(enabled_sources, str):
        enabled_sources = enabled_sources.split(",")
    normalized = {
        GENE_SOURCE_ALIASES.get(str(source).strip().lower(), str(source).strip().lower())
        for source in enabled_sources
        if str(source).strip()
    }
    unknown = normalized.difference(GENE_SOURCE_FACTORIES)
    if unknown:
        raise ValueError(f"Unknown enabled gene annotation sources: {', '.join(sorted(unknown))}")
    return normalized


def write_gene_qc(path, mag, frames):
    records = []
    for frame in frames:
        record = frame.attrs.get("annotation_qc") if frame is not None else None
        if record:
            records.append({"mag": str(mag), "level": "gene", **record})
    payload = {
        "schema_version": GENE_QC_SCHEMA,
        "mag": str(mag),
        "level": "gene",
        "sources": records,
    }
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def merge_annotations(
    gff_file,
    kegg_file,
    keggdb_file,
    keggcutoffs_file,
    pfam_file,
    ec_file,
    cazy_file,
    vf_file,
    vfdb_file,
    amr_file,
    signalp_file,
    output_file,
    card_file=None,
    defense_file=None,
    foldseek_file=None,
    foldseekdb_file=None,
    evalue_threshold=DEFAULT_EVALUE_THRESHOLD,
    identity_threshold=DEFAULT_IDENTITY_THRESHOLD,
    mag=None,
    query_coverage_threshold=DEFAULT_QUERY_COVERAGE_THRESHOLD,
    target_coverage_threshold=DEFAULT_TARGET_COVERAGE_THRESHOLD,
    enabled_sources=None,
    qc_output=None,
    ncbifam_file=None,
    ncbifam_metadata_file=None,
):
    if not mag:
        raise ValueError("MAG identity is required for the gene annotation table")

    selected = normalize_enabled_sources(enabled_sources)
    genes, gene_calls = parse_gene_calls(gff_file)
    frames = [gene_calls]
    if "kegg" in selected:
        frames.append(parse_kegg(kegg_file, keggdb_file, keggcutoffs_file, evalue_threshold))
    if "pfam" in selected:
        frames.append(parse_pfam(pfam_file, ec_file))
    if "ncbifam" in selected:
        frames.append(parse_ncbifam(ncbifam_file, ncbifam_metadata_file))
    if "cazy" in selected:
        frames.append(parse_cazy(cazy_file))
    if "virulence" in selected:
        frames.append(
            parse_vfdb(
                vf_file,
                vfdb_file,
                evalue_threshold,
                identity_threshold,
                query_coverage_threshold,
                target_coverage_threshold,
            )
        )
    if "amr" in selected:
        frames.append(parse_amr(amr_file))
    if "card" in selected:
        frames.append(parse_card(card_file))
    if "signalp" in selected:
        frames.append(parse_signalp(signalp_file))
    if "defense" in selected:
        frames.append(parse_defensefinder(defense_file))
    if "structure" in selected:
        frames.append(parse_foldseek(foldseek_file, foldseekdb_file, evalue_threshold))
    nonempty = [frame for frame in frames if frame is not None and not frame.empty]
    hits = pd.concat(nonempty, ignore_index=True) if nonempty else empty_hits()
    validate_hit_gene_ids(hits, genes, mag)

    if genes.empty:
        for column in ["contig", "start", "end", "strand"]:
            hits[column] = pd.NA
    else:
        hits = pd.merge(hits, genes, on="gene", how="left", validate="many_to_one")

    hits["mag"] = str(mag)
    for column in OUTPUT_COLUMNS:
        if column not in hits.columns:
            hits[column] = pd.NA
    hits["_source_order"] = hits["source"].map(SOURCE_ORDER).fillna(999)
    hits = hits.sort_values(["gene", "_source_order", "hit_rank"], kind="stable")
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    hits[OUTPUT_COLUMNS].to_csv(output_file, sep="\t", index=False, na_rep="")
    if qc_output:
        write_gene_qc(qc_output, mag, frames)
    return hits[OUTPUT_COLUMNS].reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description="Merge gene-level annotation hits into a lossless long-form table.")
    parser.add_argument("-gff", required=True, type=str, help="Path to the GFF file")
    parser.add_argument("-mag", "--mag", required=True, help="MAG identifier written to every output row")
    parser.add_argument("-kegg", required=False, type=str, help="Path to the KEGG HMMER table")
    parser.add_argument("-keggdb", required=False, type=str, help="Path to the KEGG hierarchy JSON")
    parser.add_argument("-keggcutoffs", required=False, type=str, help="Path to the KOfam ko_list with per-KO thresholds")
    parser.add_argument("-pfam", required=False, type=str, help="Path to the PFAM HMMER table")
    parser.add_argument("-ec", required=False, type=str, help="Path to the PFAM-to-EC mapping table")
    parser.add_argument("-ncbifam", required=False, type=str, help="Path to the NCBIfam HMMER domtblout")
    parser.add_argument(
        "-ncbifam-metadata",
        required=False,
        type=str,
        help="Path to the version-matched NCBIfam hmm_PGAP.tsv metadata",
    )
    parser.add_argument("-cazy", required=False, type=str, help="Path to dbCAN's coverage-filtered HMM result table")
    parser.add_argument("-vf", required=False, type=str, help="Path to the VFDB alignment table")
    parser.add_argument("-vfdb", required=False, type=str, help="Path to the VFDB mapping table")
    parser.add_argument("-amr", required=False, type=str, help="Path to the AMRFinderPlus report")
    parser.add_argument("-card", required=False, type=str, help="Path to the CARD/RGI protein-mode report")
    parser.add_argument("-signalp", required=False, type=str, help="Path to the SignalP table")
    parser.add_argument("-o", required=True, type=str, help="Path to the output TSV file")
    parser.add_argument("-defense", required=False, type=str, help="Path to DefenseFinder gene-level TSV")
    parser.add_argument("-foldseek", required=False, type=str, help="Path to the Foldseek/ProstT5 search output (m8)")
    parser.add_argument("-foldseekdb", required=False, type=str, help="Path to the UniProt accession -> function TSV")
    parser.add_argument(
        "--sources",
        default=None,
        help="Comma-separated enabled sources; existing files from other sources are ignored.",
    )
    parser.add_argument("--qc-output", help="Path to write per-MAG annotation QC JSON")
    parser.add_argument(
        "-evalue", "--evalue", type=float, default=DEFAULT_EVALUE_THRESHOLD,
        help="Maximum fallback e-value for applicable sources. CAZy uses dbCAN's native filters. Default: 1e-10.",
    )
    parser.add_argument(
        "-identity", "--identity", type=float, default=DEFAULT_IDENTITY_THRESHOLD,
        help="Minimum percent identity for annotation hits with identity values. Default: 50.",
    )
    parser.add_argument(
        "--query-coverage",
        type=float,
        default=DEFAULT_QUERY_COVERAGE_THRESHOLD,
        help="Minimum MMseqs query coverage as a fraction from 0 to 1. Default: 0.5.",
    )
    parser.add_argument(
        "--target-coverage",
        type=float,
        default=DEFAULT_TARGET_COVERAGE_THRESHOLD,
        help="Minimum MMseqs target coverage as a fraction from 0 to 1. Default: 0.5.",
    )
    args = parser.parse_args()
    if args.evalue < 0:
        parser.error("--evalue must be non-negative")
    if args.identity < 0 or args.identity > 100:
        parser.error("--identity must be between 0 and 100")
    if args.query_coverage < 0 or args.query_coverage > 1:
        parser.error("--query-coverage must be between 0 and 1")
    if args.target_coverage < 0 or args.target_coverage > 1:
        parser.error("--target-coverage must be between 0 and 1")

    merge_annotations(
        args.gff, args.kegg, args.keggdb, args.keggcutoffs, args.pfam, args.ec,
        args.cazy, args.vf, args.vfdb, args.amr, args.signalp, args.o,
        card_file=args.card,
        defense_file=args.defense,
        foldseek_file=args.foldseek,
        foldseekdb_file=args.foldseekdb,
        evalue_threshold=args.evalue,
        identity_threshold=args.identity,
        mag=args.mag,
        query_coverage_threshold=args.query_coverage,
        target_coverage_threshold=args.target_coverage,
        enabled_sources=args.sources,
        qc_output=args.qc_output,
        ncbifam_file=args.ncbifam,
        ncbifam_metadata_file=args.ncbifam_metadata,
    )


if __name__ == "__main__":
    main()
