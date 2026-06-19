import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SearchIO


DEFAULT_EVALUE_THRESHOLD = 1e-10
DEFAULT_IDENTITY_THRESHOLD = 50.0

# Function columns a Foldseek structural hit is allowed to fill when sequence
# homology left them empty.
STRUCTURE_FILL_COLUMNS = ["kegg", "ec", "pfam"]
# Sequence columns that, if populated, mark a gene as sequence-annotated.
SEQUENCE_EVIDENCE_COLUMNS = ["kegg", "ec", "pfam", "cazy"]


def select_lowest_evalue(group):
    return group.sort_values(by="evalue").head(1)


def select_highest_confidence(group):
    return group.sort_values(by="confidence", ascending=False).head(1)


def append_suffix_to_seqid(row):
    id_part = row["attributes"].split(";")[0].replace("ID=", "")
    suffix = id_part.split("_")[-1]
    return f"{row['seqid']}_{suffix}"


def has_content(path):
    if not path:
        return False
    path_obj = Path(path)
    return path_obj.is_file() and path_obj.stat().st_size > 0


def parse_hmmer3_tab(path):
    fields = ["accession", "bitscore", "evalue", "id", "overlap_num", "region_num"]
    # bitscore is the full-sequence score; bitscore_domain is the best single
    # domain score, needed for KOfam KOs whose curated cutoff is domain-based.
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
                domain_score = hit.hsps[0].bitscore if hit.hsps else None
                hits["bitscore_domain"].append(domain_score)

    if not query_ids:
        return pd.DataFrame(columns=columns)

    data = pd.DataFrame.from_dict(hits)
    data["gene"] = query_ids
    return data


def load_kofam_thresholds(kolist_file):
    # KOfam ships per-KO adaptive bit-score cutoffs in its ko_list file rather
    # than embedding GA/TC lines in the HMMs, so they cannot be applied with
    # hmmscan --cut_ga; they are applied here instead.
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

    kegg_rows = []
    for main in kegg_json.get("children", []):
        for broad in main.get("children", []):
            for sub in broad.get("children", []):
                for gene_node in sub.get("children", []):
                    name = gene_node.get("name", "")
                    kegg = name.split(" ")[0] if name else ""
                    description_and_ec = " ".join(name.split(" ")[1:])
                    description_parts = description_and_ec.split(" [")
                    ec = description_parts[1][:-1] if len(description_parts) > 1 else ""
                    kegg_rows.append((kegg, ec))

    return pd.DataFrame(kegg_rows, columns=["kegg", "ec"])


def parse_kegg(kegg_file, keggdb_file, kolist_file, evalue_threshold):
    kegg_df = pd.DataFrame(columns=["gene", "kegg", "ec"])
    hits = parse_hmmer3_tab(kegg_file)
    if hits.empty:
        return kegg_df

    hits["evalue"] = pd.to_numeric(hits["evalue"], errors="coerce")
    hits["bitscore"] = pd.to_numeric(hits["bitscore"], errors="coerce")
    hits["bitscore_domain"] = pd.to_numeric(hits["bitscore_domain"], errors="coerce")
    hits = hits.rename(columns={"id": "kegg"})

    thresholds = load_kofam_thresholds(kolist_file)
    if not thresholds.empty:
        # Apply each KO's curated bit-score cutoff: full-sequence score unless
        # the KO is scored on its best domain. KOs without a curated threshold
        # fall back to the flat e-value cutoff.
        hits = pd.merge(hits, thresholds, on="kegg", how="left")
        score = hits["bitscore"].where(hits["score_type"] != "domain", hits["bitscore_domain"])
        has_threshold = hits["threshold"].notna()
        passes_native = has_threshold & (score >= hits["threshold"])
        passes_fallback = ~has_threshold & (hits["evalue"] <= evalue_threshold)
        hits = hits[passes_native | passes_fallback]
    else:
        hits = hits[hits["evalue"] <= evalue_threshold]
    if hits.empty:
        return kegg_df

    hierarchy = load_kegg_hierarchy(keggdb_file)
    if not hierarchy.empty:
        hits = pd.merge(hits, hierarchy, on="kegg", how="left")
    else:
        hits["ec"] = pd.NA

    kegg_df = hits.groupby("gene", group_keys=False)[["gene", "kegg", "ec", "evalue"]].apply(
        select_lowest_evalue, include_groups=False
    ).reset_index(drop=True)
    kegg_df["ec"] = kegg_df["ec"].astype("string").str.replace("EC:", "", regex=False)
    return kegg_df


def parse_pfam(pfam_file, ec_file):
    # Pfam hits are filtered upstream by hmmscan's per-family curated gathering
    # thresholds (--cut_ga), so no flat e-value cutoff is applied here. The
    # e-value is retained only to pick the single best hit per gene.
    pfam_df = pd.DataFrame(columns=["gene", "pfam", "ec"])
    hits = parse_hmmer3_tab(pfam_file)
    if hits.empty:
        return pfam_df

    hits["evalue"] = pd.to_numeric(hits["evalue"], errors="coerce")

    hits = hits.rename(columns={"accession": "pfam"})
    hits["pfam"] = hits["pfam"].astype("string").str.split(".").str[0]
    hits = hits.groupby("gene", group_keys=False)[["gene", "pfam", "evalue"]].apply(
        select_lowest_evalue, include_groups=False
    ).reset_index(drop=True)

    if has_content(ec_file):
        pfam_to_ec = pd.read_csv(ec_file, sep="\t", comment="#", header=0)
        pfam_to_ec = pfam_to_ec[pfam_to_ec["Type"] == "GOLD"]
        pfam_to_ec = pfam_to_ec.rename(columns={
            "Confidence-Score": "confidence",
            "Pfam-Domain": "pfam",
            "EC-Number": "ec",
        })
        pfam_to_ec["confidence"] = pd.to_numeric(pfam_to_ec["confidence"], errors="coerce")
        pfam_to_ec = pfam_to_ec.groupby("pfam", group_keys=False)[["pfam", "ec", "confidence"]].apply(
            select_highest_confidence, include_groups=False
        )
        hits = pd.merge(hits, pfam_to_ec[["pfam", "ec"]], on="pfam", how="left")
    else:
        hits["ec"] = pd.NA

    return hits


def parse_cazy(cazy_file, evalue_threshold):
    cazy_df = pd.DataFrame(columns=["gene", "cazy"])
    hits = parse_hmmer3_tab(cazy_file)
    if hits.empty:
        return cazy_df

    hits["evalue"] = pd.to_numeric(hits["evalue"], errors="coerce")
    hits = hits[hits["evalue"] <= evalue_threshold]
    if hits.empty:
        return cazy_df

    hits["id"] = hits["id"].astype("string").str.replace(".hmm", "", regex=False)
    hits = hits.rename(columns={"id": "cazy"})
    cazy_df = hits.groupby("gene", group_keys=False)[["gene", "cazy", "evalue"]].apply(
        select_lowest_evalue, include_groups=False
    ).reset_index(drop=True)
    return cazy_df


def parse_amr(amr_file, amrdb_file):
    # AMR hits are filtered upstream by hmmscan's per-model trusted cutoffs
    # (--cut_tc), the same cutoffs AMRFinderPlus uses, so no flat e-value cutoff
    # is applied here. The e-value is retained only to pick the best hit per gene.
    amr_df = pd.DataFrame(columns=["gene", "resistance_type", "resistance_target"])
    hits = parse_hmmer3_tab(amr_file)
    if hits.empty:
        return amr_df

    hits["evalue"] = pd.to_numeric(hits["evalue"], errors="coerce")

    hits = hits.rename(columns={"id": "amr"})
    hits = hits.groupby("gene", group_keys=False)[["gene", "amr", "accession", "evalue"]].apply(
        select_lowest_evalue, include_groups=False
    ).reset_index(drop=True)

    if has_content(amrdb_file):
        amr_to_class = pd.read_csv(amrdb_file, sep="\t", header=0)
        amr_to_class = amr_to_class.rename(columns={"#hmm_accession": "accession"})
        hits = pd.merge(hits, amr_to_class[["accession", "subtype", "subclass"]], on="accession", how="left")
        hits = hits.rename(columns={"subtype": "resistance_type", "subclass": "resistance_target"})
    else:
        hits["resistance_type"] = pd.NA
        hits["resistance_target"] = pd.NA

    return hits[["gene", "resistance_type", "resistance_target"]]


def parse_vfdb(vf_file, vfdb_file, evalue_threshold, identity_threshold):
    vf_df = pd.DataFrame(columns=["gene", "vf", "vf_type"])
    if not has_content(vf_file):
        return vf_df

    hits = pd.read_csv(
        vf_file,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "gene", "entry", "identity", "length", "mismatches", "gaps",
            "query_start", "query_end", "target_start", "target_end", "evalue", "bitscore",
        ],
    )
    if hits.empty:
        return vf_df

    hits["evalue"] = pd.to_numeric(hits["evalue"], errors="coerce")
    hits["identity"] = pd.to_numeric(hits["identity"], errors="coerce")
    hits = hits[(hits["evalue"] <= evalue_threshold) & (hits["identity"] >= identity_threshold)]
    if hits.empty:
        return vf_df

    hits = hits.groupby("gene", group_keys=False)[["gene", "entry", "evalue"]].apply(
        select_lowest_evalue, include_groups=False
    ).reset_index(drop=True)

    if has_content(vfdb_file):
        entry_to_vf = pd.read_csv(vfdb_file, sep="\t", comment="#", header=0)
        hits = pd.merge(hits, entry_to_vf[["entry", "vf", "vf_type"]], on="entry", how="left")
    else:
        hits["vf"] = pd.NA
        hits["vf_type"] = pd.NA

    return hits[["gene", "vf", "vf_type"]]


def parse_signalp(signalp_file):
    signalp_df = pd.DataFrame(columns=["gene", "signalp"])
    if not has_content(signalp_file):
        return signalp_df

    hits = pd.read_csv(signalp_file, sep="\t", comment="#", header=None, names=["gene", "signalp", "confidence"])
    if hits.empty:
        return signalp_df

    hits["confidence"] = pd.to_numeric(hits["confidence"], errors="coerce")
    signalp_df = hits.groupby("gene", group_keys=False)[["gene", "signalp", "confidence"]].apply(
        select_highest_confidence
    ).reset_index(drop=True)
    return signalp_df[["gene", "signalp"]]


def parse_defensefinder(defense_file):
    defense_hits = pd.DataFrame(columns=["gene", "defense", "defense_type"])
    antidefense_hits = pd.DataFrame(columns=["gene", "antidefense", "antidefense_type"])
    if not has_content(defense_file):
        return defense_hits, antidefense_hits

    defense_df = pd.read_csv(defense_file, sep="\t")
    if defense_df.empty:
        return defense_hits, antidefense_hits

    defense_df = defense_df.rename(columns={"hit_id": "gene"})
    defense_df["activity"] = defense_df["activity"].fillna("")
    defense_df["gene_name"] = defense_df["gene_name"].fillna("")
    defense_df["type"] = defense_df["type"].fillna("")

    defense_hits = defense_df[defense_df["activity"] == "Defense"][["gene", "gene_name", "type"]].rename(
        columns={"gene_name": "defense", "type": "defense_type"}
    )
    antidefense_hits = defense_df[defense_df["activity"] == "Antidefense"][["gene", "gene_name", "type"]].rename(
        columns={"gene_name": "antidefense", "type": "antidefense_type"}
    )

    defense_hits = defense_hits.groupby("gene", group_keys=False).head(1).reset_index(drop=True)
    antidefense_hits = antidefense_hits.groupby("gene", group_keys=False).head(1).reset_index(drop=True)
    return defense_hits, antidefense_hits


def uniprot_accession_from_target(target):
    # AlphaFold/Swiss-Prot Foldseek targets look like "AF-P12345-F1-model_v4.cif.gz".
    match = re.match(r"AF-([^-]+)-F\d+", str(target))
    if match:
        return match.group(1)
    return str(target)


def parse_foldseek(foldseek_file, mapdb_file, evalue_threshold):
    foldseek_df = pd.DataFrame(columns=["gene", "kegg", "ec", "pfam"])
    if not has_content(foldseek_file):
        return foldseek_df

    hits = pd.read_csv(
        foldseek_file,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "gene", "target", "fident", "alnlen", "mismatch", "gapopen",
            "qstart", "qend", "tstart", "tend", "evalue", "bits",
        ],
    )
    if hits.empty:
        return foldseek_df

    hits["evalue"] = pd.to_numeric(hits["evalue"], errors="coerce")
    hits = hits[hits["evalue"] <= evalue_threshold]
    if hits.empty:
        return foldseek_df

    hits = hits.groupby("gene", group_keys=False)[["gene", "target", "evalue"]].apply(
        select_lowest_evalue, include_groups=False
    ).reset_index(drop=True)
    hits["accession"] = hits["target"].map(uniprot_accession_from_target)

    if has_content(mapdb_file):
        mapping = pd.read_csv(mapdb_file, sep="\t", comment="#", header=0)
        keep = ["accession"] + [c for c in STRUCTURE_FILL_COLUMNS if c in mapping.columns]
        hits = pd.merge(hits, mapping[keep], on="accession", how="left")

    for column in STRUCTURE_FILL_COLUMNS:
        if column not in hits.columns:
            hits[column] = pd.NA

    return hits[["gene", *STRUCTURE_FILL_COLUMNS]]


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
    amrdb_file,
    signalp_file,
    output_file,
    defense_file=None,
    foldseek_file=None,
    foldseekdb_file=None,
    evalue_threshold=DEFAULT_EVALUE_THRESHOLD,
    identity_threshold=DEFAULT_IDENTITY_THRESHOLD,
):
    annotations = pd.read_csv(
        gff_file,
        sep="\t",
        comment="#",
        header=None,
        names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"],
    )
    if annotations.empty:
        annotations = pd.DataFrame(columns=["gene", "start", "end", "strand"])
    else:
        annotations["seqid"] = annotations.apply(append_suffix_to_seqid, axis=1)
        annotations = annotations.drop(columns=["attributes", "source", "score", "type", "phase"])
        annotations = annotations.rename(columns={"seqid": "gene"})

    kegg_df = parse_kegg(kegg_file, keggdb_file, keggcutoffs_file, evalue_threshold)
    pfam_df = parse_pfam(pfam_file, ec_file)
    cazy_df = parse_cazy(cazy_file, evalue_threshold)
    amr_df = parse_amr(amr_file, amrdb_file)
    vf_df = parse_vfdb(vf_file, vfdb_file, evalue_threshold, identity_threshold)
    signalp_df = parse_signalp(signalp_file)

    annotations = pd.merge(annotations, kegg_df[["gene", "kegg", "ec"]], on="gene", how="left")
    annotations = pd.merge(annotations, pfam_df[["gene", "pfam", "ec"]], on="gene", how="left", suffixes=("", "_pfam"))
    if "ec_pfam" in annotations.columns:
        annotations["ec"] = annotations.apply(
            lambda row: row["ec_pfam"] if pd.isna(row["ec"]) or row["ec"] == "" else row["ec"],
            axis=1,
        )
        annotations.drop(columns=["ec_pfam"], inplace=True)
    annotations = pd.merge(annotations, cazy_df[["gene", "cazy"]], on="gene", how="left")
    annotations = pd.merge(annotations, amr_df[["gene", "resistance_type", "resistance_target"]], on="gene", how="left")
    annotations = pd.merge(annotations, vf_df[["gene", "vf", "vf_type"]], on="gene", how="left")
    annotations = pd.merge(annotations, signalp_df[["gene", "signalp"]], on="gene", how="left")

    defense_hits, antidefense_hits = parse_defensefinder(defense_file)
    if not defense_hits.empty:
        annotations = pd.merge(annotations, defense_hits, on="gene", how="left")
    if not antidefense_hits.empty:
        annotations = pd.merge(annotations, antidefense_hits, on="gene", how="left")

    for column in ["defense", "defense_type", "antidefense", "antidefense_type"]:
        if column not in annotations.columns:
            annotations[column] = pd.NA

    # Structural (Foldseek/ProstT5) fallback: fill kegg/ec/pfam only where
    # sequence homology left them empty, and record provenance per gene.
    def is_empty(series):
        return series.isna() | (series.astype("string").fillna("") == "")

    has_sequence = pd.Series(False, index=annotations.index)
    for column in SEQUENCE_EVIDENCE_COLUMNS:
        if column in annotations.columns:
            has_sequence |= ~is_empty(annotations[column])

    filled_by_structure = pd.Series(False, index=annotations.index)
    foldseek_df = parse_foldseek(foldseek_file, foldseekdb_file, evalue_threshold)
    if not foldseek_df.empty:
        annotations = pd.merge(annotations, foldseek_df, on="gene", how="left", suffixes=("", "_struct"))
        for column in STRUCTURE_FILL_COLUMNS:
            struct_column = f"{column}_struct"
            if struct_column not in annotations.columns:
                continue
            # Cast to object so filling an otherwise all-NaN (float) column with
            # string labels doesn't trigger a pandas dtype-incompatibility warning.
            annotations[column] = annotations[column].astype("object")
            empty = is_empty(annotations[column]) & ~is_empty(annotations[struct_column])
            annotations.loc[empty, column] = annotations.loc[empty, struct_column]
            filled_by_structure |= empty
            annotations.drop(columns=[struct_column], inplace=True)

    # Sequence provenance wins: a gene with any sequence annotation is labelled
    # "sequence" even if structure additionally filled a secondary column. In the
    # gated pipeline structure only runs on orphans, so "structure" marks genes
    # annotated purely by Foldseek.
    annotations["evidence"] = pd.NA
    annotations.loc[filled_by_structure, "evidence"] = "structure"
    annotations.loc[has_sequence, "evidence"] = "sequence"

    annotations.to_csv(output_file, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Merge gene-level annotation sources into one table.")
    parser.add_argument("-gff", required=True, type=str, help="Path to the GFF file")
    parser.add_argument("-kegg", required=False, type=str, help="Path to the KEGG HMMER table")
    parser.add_argument("-keggdb", required=False, type=str, help="Path to the KEGG hierarchy JSON")
    parser.add_argument("-keggcutoffs", required=False, type=str, help="Path to the KOfam ko_list with per-KO thresholds")
    parser.add_argument("-pfam", required=False, type=str, help="Path to the PFAM HMMER table")
    parser.add_argument("-ec", required=False, type=str, help="Path to the PFAM-to-EC mapping table")
    parser.add_argument("-cazy", required=False, type=str, help="Path to the CAZy HMMER table")
    parser.add_argument("-vf", required=False, type=str, help="Path to the VFDB alignment table")
    parser.add_argument("-vfdb", required=False, type=str, help="Path to the VFDB mapping table")
    parser.add_argument("-amr", required=False, type=str, help="Path to the AMR HMMER table")
    parser.add_argument("-amrdb", required=False, type=str, help="Path to the AMR mapping table")
    parser.add_argument("-signalp", required=False, type=str, help="Path to the SignalP table")
    parser.add_argument("-o", required=True, type=str, help="Path to the output TSV file")
    parser.add_argument("-defense", required=False, type=str, help="Path to DefenseFinder gene-level TSV")
    parser.add_argument("-foldseek", required=False, type=str, help="Path to the Foldseek/ProstT5 search output (m8)")
    parser.add_argument("-foldseekdb", required=False, type=str, help="Path to the UniProt accession -> function TSV")
    parser.add_argument(
        "-evalue",
        "--evalue",
        required=False,
        type=float,
        default=DEFAULT_EVALUE_THRESHOLD,
        help="Maximum e-value for annotation hits. Default: 1e-10.",
    )
    parser.add_argument(
        "-identity",
        "--identity",
        required=False,
        type=float,
        default=DEFAULT_IDENTITY_THRESHOLD,
        help="Minimum percent identity for annotation hits with identity values. Default: 50.",
    )

    args = parser.parse_args()
    if args.evalue < 0:
        parser.error("--evalue must be non-negative")
    if args.identity < 0 or args.identity > 100:
        parser.error("--identity must be between 0 and 100")

    merge_annotations(
        args.gff,
        args.kegg,
        args.keggdb,
        args.keggcutoffs,
        args.pfam,
        args.ec,
        args.cazy,
        args.vf,
        args.vfdb,
        args.amr,
        args.amrdb,
        args.signalp,
        args.o,
        args.defense,
        args.foldseek,
        args.foldseekdb,
        args.evalue,
        args.identity,
    )


if __name__ == "__main__":
    main()
