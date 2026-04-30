#!/usr/bin/env python3
import argparse
import gzip
import shutil
import sys
import tarfile
from pathlib import Path
from typing import Optional


BT2_SUFFIXES = {
    ".1.bt2",
    ".2.bt2",
    ".3.bt2",
    ".4.bt2",
    ".rev.1.bt2",
    ".rev.2.bt2",
}
BT2L_SUFFIXES = {
    ".1.bt2l",
    ".2.bt2l",
    ".3.bt2l",
    ".4.bt2l",
    ".rev.1.bt2l",
    ".rev.2.bt2l",
}
INDEX_SUFFIXES = tuple(sorted(BT2_SUFFIXES | BT2L_SUFFIXES, key=len, reverse=True))
FASTA_SUFFIXES = (".fna", ".fa", ".fasta")


def index_suffix(filename: str) -> Optional[str]:
    for suffix in INDEX_SUFFIXES:
        if filename.endswith(suffix):
            return suffix
    return None


def fasta_stem(filename: str) -> Optional[str]:
    name = filename
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return None


def choose_index_group(members):
    groups = {}
    for member in members:
        filename = Path(member.name).name
        suffix = index_suffix(filename)
        if not suffix:
            continue
        prefix = filename[: -len(suffix)]
        groups.setdefault(prefix, {})[suffix] = member

    complete_groups = []
    for prefix, indexed_members in groups.items():
        suffixes = set(indexed_members)
        if BT2_SUFFIXES.issubset(suffixes):
            complete_groups.append((prefix, indexed_members, BT2_SUFFIXES))
        if BT2L_SUFFIXES.issubset(suffixes):
            complete_groups.append((prefix, indexed_members, BT2L_SUFFIXES))

    if not complete_groups:
        found = ", ".join(sorted(groups)) or "none"
        raise ValueError(f"No complete Bowtie2 index found in archive. Found index basenames: {found}")
    if len(complete_groups) > 1:
        basenames = ", ".join(prefix for prefix, _, _ in complete_groups)
        raise ValueError(f"Multiple complete Bowtie2 indexes found in archive: {basenames}")
    return complete_groups[0]


def choose_fasta_member(members, reference_name: str, index_prefix: str):
    candidates = []
    for member in members:
        filename = Path(member.name).name
        stem = fasta_stem(filename)
        if stem:
            candidates.append((stem, member))

    if not candidates:
        raise ValueError("No FASTA file found in reference index archive.")

    preferred_stems = {reference_name, index_prefix}
    preferred = [member for stem, member in candidates if stem in preferred_stems]
    if len(preferred) == 1:
        return preferred[0]
    if len(candidates) == 1:
        return candidates[0][1]

    names = ", ".join(Path(member.name).name for _, member in candidates)
    raise ValueError(f"Multiple FASTA files found in archive; cannot choose automatically: {names}")


def copy_tar_member(tar, member, destination: Path, decompress_gzip: bool = False) -> None:
    source = tar.extractfile(member)
    if source is None:
        raise ValueError(f"Could not read archive member: {member.name}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source:
        if decompress_gzip:
            with gzip.GzipFile(fileobj=source) as decompressed, destination.open("wb") as output:
                shutil.copyfileobj(decompressed, output)
        else:
            with destination.open("wb") as output:
                shutil.copyfileobj(source, output)


def extract_reference_index(archive: Path, reference_name: str, output_dir: Path) -> None:
    if not tarfile.is_tarfile(str(archive)):
        raise ValueError(f"Reference index is not a readable tar archive: {archive}")

    output_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(str(archive)) as tar:
        members = [member for member in tar.getmembers() if member.isfile()]
        index_prefix, index_members, required_suffixes = choose_index_group(members)
        fasta_member = choose_fasta_member(members, reference_name, index_prefix)

        fasta_name = Path(fasta_member.name).name
        copy_tar_member(
            tar,
            fasta_member,
            output_dir / f"{reference_name}.fna",
            decompress_gzip=fasta_name.endswith(".gz"),
        )

        for suffix in sorted(required_suffixes):
            copy_tar_member(tar, index_members[suffix], output_dir / f"{reference_name}{suffix}")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Extract a reference FASTA and Bowtie2 index files from a tar archive.")
    parser.add_argument("--archive", required=True, type=Path, help="Reference index tarball")
    parser.add_argument("--reference", required=True, help="Output reference basename")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for extracted reference files")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        extract_reference_index(args.archive, args.reference, args.output_dir)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
