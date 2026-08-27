"""Resolve assembly inputs for the :command:`drakkar amr` workflow.

The AMR workflow operates on assemblies rather than MAGs.  Its manifest keeps
the user-facing assembly identity separate from the FASTA basename and carries
the small amount of metadata that changes how the callers are run.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
from pathlib import Path

from drakkar.input_errors import InputFileError, report_input_resolution_errors
from drakkar.input_tables import read_input_table


FASTA_RE = re.compile(r"\.(?:fa|fna|fasta)(?:\.gz)?$", re.IGNORECASE)
ASSEMBLY_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
ASSEMBLY_ID_COLUMNS = ("assembly_id", "assembly", "sample")
ASSEMBLY_PATH_COLUMNS = ("assembly_path", "path", "fasta")
ASSEMBLY_TYPES = {"metagenome", "isolate"}
IUPAC_DNA = set("ACGTRYSWKMBDHVN")


def _first_column(columns, candidates):
    return next((candidate for candidate in candidates if candidate in columns), None)


def _table_value(value):
    """Return a stripped table cell, treating pandas' NaN as missing."""
    if value is None:
        return ""
    try:
        if value != value:  # NaN is the only ordinary scalar unequal to itself.
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _open_fasta(path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inspect_assembly(path):
    """Validate a nucleotide FASTA and return stable input statistics."""
    path = Path(path).resolve()
    if not path.is_file():
        raise InputFileError(f"Assembly FASTA not found: {path}")
    if path.stat().st_size == 0:
        raise InputFileError(f"Assembly FASTA is empty: {path}")
    if not FASTA_RE.search(path.name):
        raise InputFileError(
            f"Assembly file must end in .fa, .fna, .fasta, or the corresponding .gz form: {path}"
        )

    identifiers = set()
    current = None
    contig_count = 0
    total_length = 0
    try:
        with _open_fasta(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if line.startswith(">"):
                    identifier = line[1:].strip().split(None, 1)[0]
                    if not identifier:
                        raise InputFileError(
                            f"Assembly FASTA has an empty identifier at {path}:{line_number}."
                        )
                    if identifier in identifiers:
                        raise InputFileError(
                            f"Assembly FASTA contains duplicate first-token contig identifier "
                            f"{identifier!r}: {path}"
                        )
                    identifiers.add(identifier)
                    current = identifier
                    contig_count += 1
                    continue
                sequence = "".join(line.split())
                if not sequence:
                    continue
                if current is None:
                    raise InputFileError(
                        f"Assembly FASTA contains sequence before its first header at {path}:{line_number}."
                    )
                invalid = sorted(set(sequence.upper()) - IUPAC_DNA)
                if invalid:
                    preview = "".join(invalid[:10])
                    raise InputFileError(
                        f"Assembly FASTA contains non-IUPAC nucleotide character(s) "
                        f"{preview!r} at {path}:{line_number}."
                    )
                total_length += len(sequence)
    except (OSError, EOFError) as exc:
        raise InputFileError(f"Could not read assembly FASTA {path}: {exc}") from exc

    if contig_count == 0 or total_length == 0:
        raise InputFileError(f"Assembly FASTA contains no nucleotide sequences: {path}")
    return {
        "path": str(path),
        "sha256": _sha256(path),
        "size_bytes": path.stat().st_size,
        "contig_count": contig_count,
        "total_length": total_length,
    }


def _normalize_assembly_id(value, *, context):
    assembly_id = _table_value(value)
    if not assembly_id:
        raise InputFileError(f"Missing assembly_id {context}.")
    if not ASSEMBLY_ID_RE.fullmatch(assembly_id):
        raise InputFileError(
            f"Invalid assembly_id {assembly_id!r} {context}; use letters, numbers, '.', '_', and '-' only."
        )
    return assembly_id


def _normalize_type(value, default, *, context):
    assembly_type = (_table_value(value) or default).strip().lower()
    if assembly_type not in ASSEMBLY_TYPES:
        raise InputFileError(
            f"Invalid assembly_type {assembly_type!r} {context}; expected metagenome or isolate."
        )
    return assembly_type


def _record(assembly_id, path, assembly_type, organism="", source_layout="manifest"):
    stats = inspect_assembly(path)
    return {
        **stats,
        "assembly_type": assembly_type,
        "organism": _table_value(organism),
        "source_layout": source_layout,
    }


def assemblies_from_table(table_path, default_type="metagenome"):
    """Read a CSV/TSV assembly manifest into the workflow mapping."""
    table_path = Path(table_path).resolve()
    frame = read_input_table(
        str(table_path), label="assembly manifest", dtype=str
    )
    id_column = _first_column(frame.columns, ASSEMBLY_ID_COLUMNS)
    path_column = _first_column(frame.columns, ASSEMBLY_PATH_COLUMNS)
    missing = []
    if id_column is None:
        missing.append("assembly_id")
    if path_column is None:
        missing.append("assembly_path")
    if missing:
        raise InputFileError(
            f"Assembly manifest {table_path} is missing required column(s): {', '.join(missing)}."
        )

    assemblies = {}
    errors = []
    for index, row in frame.iterrows():
        context = f"on row {index + 1} of {table_path}"
        try:
            assembly_id = _normalize_assembly_id(row.get(id_column), context=context)
            if assembly_id in assemblies:
                raise InputFileError(f"Duplicate assembly_id {assembly_id!r} {context}.")
            raw_path = _table_value(row.get(path_column))
            if not raw_path:
                raise InputFileError(f"Missing assembly_path {context}.")
            path = Path(raw_path).expanduser()
            if not path.is_absolute():
                path = table_path.parent / path
            assembly_type = _normalize_type(
                row.get("assembly_type"), default_type, context=context
            )
            assemblies[assembly_id] = _record(
                assembly_id,
                path,
                assembly_type,
                row.get("organism", ""),
            )
        except InputFileError as exc:
            errors.append(str(exc))

    if not assemblies and not errors:
        errors.append(f"No assembly rows were found in {table_path}.")
    if errors:
        raise InputFileError("\n".join(errors))
    return assemblies


def _fasta_files(directory):
    return sorted(
        path for path in directory.iterdir()
        if path.is_file() and FASTA_RE.search(path.name)
    )


def _choose_nested_fasta(directory):
    preferred = [
        directory / f"{directory.name}.fna",
        directory / f"{directory.name}.fna.gz",
        directory / "final.contigs.fa",
        directory / "final.contigs.fa.gz",
    ]
    for candidate in preferred:
        if candidate.is_file():
            return candidate
    candidates = _fasta_files(directory)
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        names = ", ".join(path.name for path in candidates)
        raise InputFileError(
            f"Assembly directory {directory} contains multiple FASTA files and no unambiguous "
            f"<directory>.fna or final.contigs.fa candidate: {names}"
        )
    return None


def assemblies_from_directory(input_dir, default_type="metagenome"):
    """Discover flat FASTAs or one assembly FASTA per immediate child folder."""
    root = Path(input_dir).resolve()
    if not root.is_dir():
        raise InputFileError(f"Assembly input directory not found: {root}")

    drakkar_megahit = root / "cataloging" / "megahit"
    if drakkar_megahit.is_dir():
        root = drakkar_megahit
        source_layout = "drakkar_cataloging_megahit"
    else:
        source_layout = "directory"

    candidates = []
    candidates.extend(
        (FASTA_RE.sub("", path.name), path, source_layout)
        for path in _fasta_files(root)
    )
    for child in sorted(path for path in root.iterdir() if path.is_dir()):
        selected = _choose_nested_fasta(child)
        if selected is not None:
            layout = (
                source_layout
                if source_layout == "drakkar_cataloging_megahit"
                else "nested_directory"
            )
            candidates.append((child.name, selected, layout))

    if not candidates:
        raise InputFileError(
            f"No assembly FASTA files were found in {root}. Expected a flat FASTA directory or "
            "one <assembly>/<assembly>.fna or <assembly>/final.contigs.fa file per child folder."
        )

    assemblies = {}
    errors = []
    for raw_id, path, layout in candidates:
        try:
            assembly_id = _normalize_assembly_id(raw_id, context=f"derived from {path}")
            if assembly_id in assemblies:
                raise InputFileError(
                    f"Duplicate assembly_id {assembly_id!r} discovered in {root}."
                )
            assemblies[assembly_id] = _record(
                assembly_id, path, default_type, source_layout=layout
            )
        except InputFileError as exc:
            errors.append(str(exc))
    if errors:
        raise InputFileError("\n".join(errors))
    return assemblies


def write_amr_assemblies_manifest(*, output, input_dir=None, table=None, default_type="metagenome"):
    """Resolve one input mode and write ``data/amr_assemblies.json``."""
    try:
        if table:
            assemblies = assemblies_from_table(table, default_type=default_type)
        elif input_dir:
            assemblies = assemblies_from_directory(input_dir, default_type=default_type)
        else:
            raise InputFileError("Provide either an assembly directory (-i) or manifest (-f).")
    except InputFileError as exc:
        report_input_resolution_errors(str(exc).splitlines())

    output_path = Path(output) / "data" / "amr_assemblies.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(assemblies, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return assemblies
