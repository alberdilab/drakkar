#!/usr/bin/env python3
"""Validate one versioned NCBIfam/PGAP profile release before hmmpress."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


REQUIRED_METADATA_COLUMNS = {
    "ncbi_accession",
    "label",
    "sequence_cutoff",
    "domain_cutoff",
    "family_type",
    "source",
}


def read_profile_cutoffs(metadata_path):
    """Return metadata by exact accession, rejecting missing native cutoffs."""
    metadata_path = Path(metadata_path)
    with metadata_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames:
            reader.fieldnames = [name.lstrip("#") for name in reader.fieldnames]
        missing = REQUIRED_METADATA_COLUMNS.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "NCBIfam metadata is missing required column(s): "
                + ", ".join(sorted(missing))
            )

        records = {}
        for row in reader:
            accession = str(row.get("ncbi_accession") or "").strip()
            # hmm_PGAP.tsv also describes Pfam models, but NCBI does not ship
            # those HMMs in hmm_PGAP.LIB. Only exact NCBI/TIGR accessions can
            # correspond to an installed profile.
            if not re.fullmatch(r"(?:NF\d+|TIGR\d+)\.\d+", accession):
                continue
            if accession in records:
                raise ValueError(f"Duplicate NCBIfam metadata accession: {accession}")
            if not str(row.get("sequence_cutoff") or "").strip() or not str(
                row.get("domain_cutoff") or ""
            ).strip():
                raise ValueError(
                    f"NCBIfam profile {accession} has no complete trusted cutoff; "
                    "Drakkar does not use an E-value fallback."
                )
            records[accession] = row
    if not records:
        raise ValueError("NCBIfam metadata contains no versioned NF or TIGR profiles")
    return records


def read_library_accessions(library_path):
    """Return library accessions and fail if a profile omits ACC or TC."""
    accessions = set()
    current_accession = None
    has_tc = False
    model_count = 0

    def finish_model():
        if not model_count:
            return
        if not current_accession:
            raise ValueError(f"NCBIfam profile {model_count} has no ACC field")
        if not has_tc:
            raise ValueError(
                f"NCBIfam profile {current_accession} has no trusted cutoff; "
                "Drakkar does not use an E-value fallback."
            )
        if current_accession in accessions:
            raise ValueError(f"Duplicate NCBIfam library accession: {current_accession}")
        accessions.add(current_accession)

    with Path(library_path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("HMMER"):
                finish_model()
                model_count += 1
                current_accession = None
                has_tc = False
            elif line.startswith("ACC"):
                parts = line.split()
                current_accession = parts[1] if len(parts) > 1 else None
            elif line.startswith("TC"):
                values = line.split()[1:3]
                has_tc = len(values) == 2 and all(value not in {"", "-"} for value in values)
            elif line.startswith("//"):
                finish_model()
                model_count = 0
                current_accession = None
                has_tc = False

    finish_model()
    if not accessions:
        raise ValueError("NCBIfam profile library contains no models")
    return accessions


def validate_release(library, metadata, release_notes, expected_release):
    notes = Path(release_notes).read_text(encoding="utf-8")
    match = re.search(r"Release number/name:\s*hmm_PGAP/(\d+\.\d+)", notes)
    actual_release = match.group(1) if match else None
    if actual_release != expected_release:
        raise ValueError(
            f"NCBIfam release notes identify {actual_release or 'no release'}, "
            f"expected {expected_release}"
        )

    metadata_records = read_profile_cutoffs(metadata)
    library_accessions = read_library_accessions(library)
    missing_metadata = sorted(library_accessions.difference(metadata_records))
    if missing_metadata:
        preview = ", ".join(missing_metadata[:10])
        raise ValueError(
            f"NCBIfam metadata has no row for {len(missing_metadata)} installed profile(s): {preview}"
        )
    return {
        "release": actual_release,
        "profiles": len(library_accessions),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--release-notes", required=True)
    parser.add_argument("--expect-release", required=True)
    args = parser.parse_args()
    result = validate_release(
        args.library, args.metadata, args.release_notes, args.expect_release
    )
    print(
        f"Validated NCBIfam/PGAP release {result['release']} "
        f"with {result['profiles']} trusted-cutoff profiles."
    )


if __name__ == "__main__":
    main()
