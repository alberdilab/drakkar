"""Section definitions and the capability probe for the Drakkar report.

Drakkar is modular, so a missing source file is the normal case rather than an
error. The probe reports which sections a given output directory can actually
support, and the report command renders the intersection of what is available
and what the user asked for.
"""

from pathlib import Path

# Section -> ordered source files. A section is available when at least one of
# its `required` files exists and is non-empty; `optional` files enrich it.
SECTION_SOURCES = {
    "preprocessing": {
        "label": "Preprocessing",
        "required": ["preprocessing.tsv"],
        "optional": [],
    },
    "cataloging": {
        "label": "Cataloging",
        "required": ["cataloging.tsv"],
        "optional": [],
    },
    "dereplication": {
        "label": "Dereplication",
        "required": ["dereplicating.tsv"],
        # dRep's own data tables: the cluster assignments, the within-cluster
        # ANI comparisons and the raw MASH distances. They turn the two-number
        # summary into a picture of how the ANI threshold acted.
        "optional": [
            "dereplicating/drep/data_tables/Cdb.csv",
            "dereplicating/drep/data_tables/Wdb.csv",
            "dereplicating/drep/data_tables/Ndb.csv",
            "dereplicating/drep/data_tables/Mdb.csv",
        ],
    },
    "profiling": {
        "label": "Profiling",
        "required": [
            "profiling_genomes/final/counts.tsv",
            "profiling_genomes/final/mags.tsv",
        ],
        "optional": ["profiling_genomes/final/bases.tsv"],
    },
    "taxonomy": {
        "label": "Taxonomy",
        "required": ["annotating/genome_taxonomy.tsv"],
        "optional": [],
    },
    "function": {
        "label": "Functional annotation",
        "required": [
            "annotating/gene_annotations.tsv.xz",
            "annotating/cluster_annotations.tsv.xz",
        ],
        "optional": ["annotating/annotation_qc.tsv"],
    },
    "expression": {
        "label": "Expression",
        "required": ["expressing/gene_counts.tsv.xz"],
        "optional": [],
    },
    # The resources section has no fixed file names: run metadata and the
    # benchmark artefacts are all stamped with a run id, so they are discovered
    # by glob in `probe_section` rather than listed here.
    "resources": {
        "label": "Runs and resources",
        "required": [],
        "optional": [],
    },
}

# Per-assembly Binette reports, one ``<assembly>.tsv`` per assembly.
BIN_REPORT_DIRECTORY = "cataloging/final"

# Written by drakkar.benchmark next to the run metadata; absent for runs that
# were not launched on SLURM.
BENCHMARK_DIRECTORY = "benchmark"
BENCHMARK_SUMMARY_SUFFIX = "_resources.yaml"

SECTION_ORDER = (
    "preprocessing",
    "cataloging",
    "dereplication",
    "profiling",
    "taxonomy",
    "function",
    "expression",
    "resources",
)

ALL_SECTIONS = frozenset(SECTION_ORDER)


class SectionError(ValueError):
    """Raised when the user requests a section that does not exist."""


def parse_sections(value):
    """Parse a comma-separated --sections value into an ordered tuple.

    ``None`` or ``all`` selects every section. Unknown names raise
    ``SectionError`` naming the valid choices, rather than being ignored.
    """
    if value is None:
        return SECTION_ORDER
    names = [item.strip().lower() for item in str(value).split(",") if item.strip()]
    if not names:
        return SECTION_ORDER
    if "all" in names:
        return SECTION_ORDER
    unknown = [name for name in names if name not in ALL_SECTIONS]
    if unknown:
        raise SectionError(
            f"Unknown report section(s): {', '.join(sorted(set(unknown)))}. "
            f"Valid sections are: {', '.join(SECTION_ORDER)}, all."
        )
    # Preserve the canonical order rather than the order the user typed.
    selected = set(names)
    return tuple(name for name in SECTION_ORDER if name in selected)


def _is_usable(path):
    """A source file counts only when it exists and holds more than a header."""
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def find_run_metadata(output_dir):
    """Return the run metadata YAML files in an output directory, oldest first.

    The benchmark roll-ups live beside them and match the same glob, so they
    are excluded here by suffix: they describe a run's resource usage, not the
    run itself, and would otherwise overwrite its provenance row.
    """
    output_path = Path(output_dir)
    try:
        candidates = sorted(output_path.glob("drakkar_*.yaml"))
    except OSError:
        return []
    return [
        path for path in candidates
        if not path.name.endswith(BENCHMARK_SUMMARY_SUFFIX)
    ]


def find_bin_reports(output_dir):
    """Return the per-assembly Binette bin quality reports, in assembly order.

    ``cataloging/final/<assembly>.tsv`` is a verbatim copy of Binette's
    ``final_bins_quality_reports.tsv``, so it carries the ``origin`` column
    that says which binner each final bin came from — or that Binette built it
    itself. The file names are assembly ids rather than a fixed name, so they
    are found by glob rather than listed in ``SECTION_SOURCES``.
    """
    directory = Path(output_dir) / BIN_REPORT_DIRECTORY
    try:
        return sorted(path for path in directory.glob("*.tsv") if _is_usable(path))
    except OSError:
        return []


def find_benchmark_summaries(output_dir):
    """Return the per-run benchmark roll-up YAMLs, oldest first."""
    output_path = Path(output_dir)
    try:
        return sorted(output_path.glob(f"drakkar_*{BENCHMARK_SUMMARY_SUFFIX}"))
    except OSError:
        return []


def find_benchmark_tables(output_dir, kind):
    """Return the per-run ``jobs`` or ``rules`` benchmark TSVs, oldest first."""
    output_path = Path(output_dir) / BENCHMARK_DIRECTORY
    try:
        return sorted(output_path.glob(f"drakkar_*.{kind}.tsv"))
    except OSError:
        return []


def benchmark_run_id(path, kind=None):
    """Recover the run id a benchmark file belongs to from its name."""
    name = Path(path).name
    if kind is not None:
        stem = name[: -len(f".{kind}.tsv")] if name.endswith(f".{kind}.tsv") else name
    elif name.endswith(BENCHMARK_SUMMARY_SUFFIX):
        stem = name[: -len(BENCHMARK_SUMMARY_SUFFIX)]
    else:
        stem = Path(name).stem
    return stem[len("drakkar_"):] if stem.startswith("drakkar_") else stem


def probe_section(output_dir, section):
    """Return availability details for one section."""
    output_path = Path(output_dir)
    spec = SECTION_SOURCES[section]

    if section == "resources":
        runs = find_run_metadata(output_path)
        benchmarks = [
            path
            for group in (
                find_benchmark_summaries(output_path),
                find_benchmark_tables(output_path, "jobs"),
                find_benchmark_tables(output_path, "rules"),
            )
            for path in group
        ]
        present = runs + benchmarks
        missing = [] if runs else ["drakkar_<run_id>.yaml"]
        if runs and not benchmarks:
            # Not an error: only SLURM runs are benchmarked, so this names what
            # the section will be missing rather than blocking it.
            missing.append(f"drakkar_<run_id>{BENCHMARK_SUMMARY_SUFFIX}")
        return {
            "section": section,
            "label": spec["label"],
            "available": bool(present),
            "present": [str(path.relative_to(output_path)) for path in present],
            "missing": missing,
        }

    present = []
    missing = []
    for relative in spec["required"]:
        candidate = output_path / relative
        (present if _is_usable(candidate) else missing).append(relative)
    optional_present = [
        relative
        for relative in spec["optional"]
        if _is_usable(output_path / relative)
    ]
    if section == "cataloging":
        optional_present.extend(
            str(path.relative_to(output_path))
            for path in find_bin_reports(output_path)
        )
    return {
        "section": section,
        "label": spec["label"],
        "available": bool(present),
        "present": present + optional_present,
        "missing": missing,
    }


def probe(output_dir, sections=None):
    """Probe an output directory for every requested section.

    Returns a list of per-section dictionaries in canonical order, each with
    ``available``, ``present`` and ``missing`` keys.
    """
    requested = sections if sections is not None else SECTION_ORDER
    return [probe_section(output_dir, section) for section in requested]


def available_sections(output_dir, sections=None):
    """Return only the section names that can actually be rendered."""
    return tuple(
        result["section"] for result in probe(output_dir, sections) if result["available"]
    )
