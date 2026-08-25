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
        "optional": [],
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
    "resources": {
        "label": "Runs and resources",
        "required": [],
        "optional": [],
    },
}

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
    """Return the run metadata YAML files in an output directory, oldest first."""
    output_path = Path(output_dir)
    try:
        return sorted(output_path.glob("drakkar_*.yaml"))
    except OSError:
        return []


def probe_section(output_dir, section):
    """Return availability details for one section."""
    output_path = Path(output_dir)
    spec = SECTION_SOURCES[section]

    if section == "resources":
        runs = find_run_metadata(output_path)
        return {
            "section": section,
            "label": spec["label"],
            "available": bool(runs),
            "present": [str(path.relative_to(output_path)) for path in runs],
            "missing": [] if runs else ["drakkar_<run_id>.yaml"],
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
