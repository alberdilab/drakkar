#!/usr/bin/env python3
"""Column vocabulary shared by the two consumers of AMRFinderPlus and RGI output.

``amr_digest.py`` turns those outputs into assembly-level loci and
``merge_gene_annotations.py`` turns them into gene-level annotation rows. Both
read the same native tables, and both have to survive the column renames the
two tools have made across releases, so the alias groups live here once rather
than in each parser.

Aliases are matched case-insensitively and ignoring punctuation, so
``"% Identity to reference"`` and ``"%Identity_to_reference"`` resolve to the
same field.
"""

from __future__ import annotations

import re


# AMRFinderPlus renamed several columns in 4.0 ("Gene symbol" became "Element
# symbol", "Element type" became "Type"). Each tuple lists the current name
# first and the names Drakkar still has to read second.
AMRFINDER_COLUMNS = {
    "gene": ("Protein identifier", "Protein id"),
    "contig": ("Contig id", "Contig"),
    "start": ("Start",),
    "end": ("Stop", "End"),
    "strand": ("Strand",),
    "symbol": ("Element symbol", "Gene symbol"),
    "name": ("Element name", "Sequence name"),
    "type": ("Type", "Element type"),
    "subtype": ("Subtype", "Element subtype"),
    "method": ("Method",),
    "drug_class": ("Class",),
    "drug_subclass": ("Subclass",),
    "identity": ("% Identity to reference", "% Identity to reference sequence"),
    "reference_coverage": ("% Coverage of reference", "% Coverage of reference sequence"),
    "alignment_length": ("Alignment length",),
    "reference_accession": ("Closest reference accession", "Accession of closest sequence"),
    "reference_name": ("Closest reference name", "Name of closest sequence"),
    "hmm_accession": ("HMM accession", "HMM id"),
    "hmm_description": ("HMM description",),
    "hierarchy_node": ("Hierarchy node",),
}

# The columns AMRFinderPlus output must carry for Drakkar to parse it at all.
AMRFINDER_REQUIRED = (
    AMRFINDER_COLUMNS["gene"],
    AMRFINDER_COLUMNS["symbol"],
    AMRFINDER_COLUMNS["type"],
    AMRFINDER_COLUMNS["method"],
)

RGI_COLUMNS = {
    "gene": ("ORF_ID", "ORF ID"),
    "contig": ("Contig",),
    "start": ("Start",),
    "end": ("Stop", "End"),
    "strand": ("Orientation", "Strand"),
    "aro_name": ("Best_Hit_ARO", "Best Hit ARO"),
    "aro": ("ARO",),
    "cutoff": ("Cut_Off", "Cut Off"),
    "model_type": ("Model_type", "Model type"),
    "drug_class": ("Drug Class",),
    "drug_subclass": ("Antibiotic",),
    "resistance_mechanism": ("Resistance Mechanism",),
    "gene_family": ("AMR Gene Family",),
    "identity": ("Best_Identities", "Best Identities"),
    "reference_coverage": ("Percentage Length of Reference Sequence",),
    "bitscore": ("Best_Hit_bit-score", "Best Hit bit-score"),
    "threshold": ("Pass_bit-score", "Pass bit-score"),
    "snps": ("SNPs_in_Best_Hit_ARO", "SNPs in Best Hit ARO"),
    "note": ("Note",),
}

RGI_REQUIRED = (
    RGI_COLUMNS["gene"],
    RGI_COLUMNS["aro_name"],
    RGI_COLUMNS["cutoff"],
)

# RGI reports its own confidence tier rather than a single cutoff. Perfect is an
# exact match to a curated reference, Strict is above the curated bitscore
# cutoff for the model, and Loose is below it and off by default.
RGI_CUTOFF_RANK = {"PERFECT": 3, "STRICT": 2, "LOOSE": 1}

# AMRFinderPlus's own evidence ranking, best first, from its Methods page:
# ALLELE > EXACT > BLAST > INTERNAL_STOP > PARTIAL_CONTIG_END > PARTIAL > HMM.
# Drakkar ranks a gene's competing AMR calls the same way so that hit_rank 1 is
# the call AMRFinderPlus itself would prefer.
AMRFINDER_METHOD_RANK = {
    "ALLELEX": 8, "ALLELEP": 8, "ALLELE": 8,
    "EXACTX": 7, "EXACTP": 7, "EXACT": 7,
    "BLASTX": 6, "BLASTP": 6, "BLAST": 6,
    "INTERNALSTOP": 5,
    "PARTIALCONTIGENDX": 4, "PARTIALCONTIGENDP": 4, "PARTIALCONTIGEND": 4,
    "PARTIALX": 3, "PARTIALP": 3, "PARTIAL": 3,
    "HMM": 2,
    "POINTX": 1, "POINTP": 1, "POINTN": 1, "POINT": 1,
}


def canonical_key(value):
    """Normalize a column name so punctuation and case differences do not matter."""
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def clean(value):
    """Return a stripped string, mapping the usual null spellings to ''."""
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() in {"", "na", "nan", "none", "null"} else value


def indexed_row(row):
    """Index one native record by canonical column name."""
    return {canonical_key(key): value for key, value in row.items()}


def field(indexed, aliases, default=""):
    """Read the first alias present and non-empty in an indexed row."""
    for alias in aliases:
        value = clean(indexed.get(canonical_key(alias)))
        if value:
            return value
    return default


def missing_columns(fieldnames, required):
    """Return the alias groups that no column in fieldnames satisfies."""
    present = {canonical_key(name) for name in fieldnames}
    return [
        "/".join(group)
        for group in required
        if not any(canonical_key(alias) in present for alias in group)
    ]


def method_rank(method):
    """Score an AMRFinderPlus method string; higher is stronger evidence."""
    return AMRFINDER_METHOD_RANK.get(canonical_key(method).upper(), 0)


def cutoff_rank(cutoff):
    """Score an RGI Cut_Off tier; higher is stronger evidence."""
    return RGI_CUTOFF_RANK.get(canonical_key(cutoff).upper(), 0)
