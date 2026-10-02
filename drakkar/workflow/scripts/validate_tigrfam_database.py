#!/usr/bin/env python3
"""Validate a legacy TIGRFAM profile release before hmmpress."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def read_library_accessions(library_path):
    """Return exact legacy accessions and reject profiles without TC1/TC2."""
    accessions = set()
    current_accession = None
    trusted_cutoffs = None
    in_model = False

    def finish_model():
        if not in_model:
            return
        if not current_accession:
            raise ValueError("TIGRFAM profile has no ACC field")
        if not re.fullmatch(r"TIGR\d{5}", current_accession):
            raise ValueError(
                f"TIGRFAM profile has a non-legacy accession: {current_accession}"
            )
        if trusted_cutoffs is None:
            raise ValueError(
                f"TIGRFAM profile {current_accession} has no complete trusted cutoff; "
                "Drakkar does not use an E-value fallback."
            )
        if current_accession in accessions:
            raise ValueError(f"Duplicate TIGRFAM library accession: {current_accession}")
        accessions.add(current_accession)

    with Path(library_path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("HMMER"):
                finish_model()
                in_model = True
                current_accession = None
                trusted_cutoffs = None
            elif line.startswith("ACC"):
                parts = line.split()
                current_accession = parts[1] if len(parts) > 1 else None
            elif line.startswith("TC"):
                values = [value.rstrip(";") for value in line.split()[1:3]]
                if len(values) == 2:
                    try:
                        trusted_cutoffs = tuple(float(value) for value in values)
                    except ValueError:
                        trusted_cutoffs = None
            elif line.startswith("//"):
                finish_model()
                in_model = False
                current_accession = None
                trusted_cutoffs = None

    finish_model()
    if not accessions:
        raise ValueError("TIGRFAM profile library contains no models")
    return accessions


def validate_release(library, release_notes, expected_release):
    notes = Path(release_notes).read_text(encoding="utf-8", errors="replace")
    if expected_release not in notes:
        raise ValueError(
            f"TIGRFAM release notes do not identify requested release {expected_release}"
        )
    accessions = read_library_accessions(library)
    return {"release": expected_release, "profiles": len(accessions)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", required=True)
    parser.add_argument("--release-notes", required=True)
    parser.add_argument("--expect-release", required=True)
    args = parser.parse_args()
    result = validate_release(args.library, args.release_notes, args.expect_release)
    print(
        f"Validated legacy TIGRFAM release {result['release']} "
        f"with {result['profiles']} trusted-cutoff profiles."
    )


if __name__ == "__main__":
    main()
