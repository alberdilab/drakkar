"""The `drakkar report` command.

Builds ``drakkar.db``, a lean SQLite projection of whatever a Drakkar output
directory happens to contain. The HTML report is rendered from that database in
a later step, so the database is always the single source of truth and the
source tables are read exactly once.
"""

from pathlib import Path

from drakkar import __version__
from drakkar.cli_context import ERROR, INFO, RESET
from drakkar.output import print, section
from drakkar.report import ingest as _ingest
from drakkar.report.schema import (
    SCHEMA_VERSION,
    connect,
    create_schema,
    finalize,
    read_schema_version,
)
from drakkar.report.sources import (
    SECTION_SOURCES,
    SectionError,
    parse_sections,
    probe,
)

DATABASE_NAME = "drakkar.db"
REPORT_NAME = "drakkar_report.html"


def build_database(output_dir, sections, db_path, primary_hits_only=False):
    """Ingest every available requested section into a fresh database.

    Returns a list of per-section result dictionaries.
    """
    connection = connect(db_path)
    try:
        create_schema(connection, drakkar_version=__version__)
        results = []
        for entry in probe(output_dir, sections):
            name = entry["section"]
            if not entry["available"]:
                results.append({**entry, "rows": None})
                continue
            loader = _ingest.SECTION_LOADERS[name]
            if name == "function":
                rows = loader(connection, output_dir, primary_hits_only=primary_hits_only)
            else:
                rows = loader(connection, output_dir)
            results.append({**entry, "rows": rows})
        finalize(connection)
        return results
    finally:
        connection.close()


def _check_existing_database(db_path, force):
    """Decide whether an existing database can be reused, rebuilt, or blocks."""
    if not db_path.exists():
        return True
    if force:
        try:
            db_path.unlink()
        except OSError as exc:
            print(f"{ERROR}ERROR:{RESET} Cannot remove existing report database: {db_path}")
            print(f"{exc.__class__.__name__}: {exc}")
            return False
        # WAL sidecars from an interrupted build would otherwise be reapplied.
        for suffix in ("-wal", "-shm"):
            sidecar = Path(str(db_path) + suffix)
            if sidecar.exists():
                try:
                    sidecar.unlink()
                except OSError:
                    pass
        return True

    connection = None
    try:
        connection = connect(db_path)
        stored = read_schema_version(connection)
    except Exception:
        stored = None
    finally:
        if connection is not None:
            connection.close()

    if stored != SCHEMA_VERSION:
        print(f"{ERROR}ERROR:{RESET} Existing report database uses schema version {stored}, "
              f"but this Drakkar build expects version {SCHEMA_VERSION}: {db_path}")
        print("Rebuild it from the source tables with 'drakkar report --force'.")
        return False
    return True


def run_report(output_dir, sections=None, db_only=False, force=False, primary_hits_only=False):
    """Entry point for `drakkar report`."""
    section("BUILDING DRAKKAR REPORT")

    output_path = Path(output_dir)
    if not output_path.is_dir():
        print(f"{ERROR}ERROR:{RESET} Output directory not found: {output_path}")
        return 1

    try:
        requested = parse_sections(sections)
    except SectionError as exc:
        print(f"{ERROR}ERROR:{RESET} {exc}")
        return 1

    db_path = output_path / DATABASE_NAME
    if not _check_existing_database(db_path, force):
        return 1

    results = build_database(
        output_path, requested, db_path, primary_hits_only=primary_hits_only
    )

    rendered = [entry for entry in results if entry["available"]]
    skipped = [entry for entry in results if not entry["available"]]

    if not rendered:
        print(f"{ERROR}ERROR:{RESET} No report sections could be built from: {output_path}")
        print("None of the expected Drakkar output tables were found there.")
        _print_skipped(skipped)
        return 1

    print(f"{INFO}INFO:{RESET} Report database written: {db_path}")
    for entry in rendered:
        rows = entry["rows"]
        detail = f"{rows} rows" if rows is not None else "no rows"
        print(f"  {entry['label']:<24} {detail}")

    # Sections the user explicitly asked for but whose inputs are absent are
    # named individually; an unfiltered run just lists what was skipped.
    _print_skipped(skipped, explicit=sections is not None)

    if not db_only:
        print(f"{INFO}INFO:{RESET} HTML rendering is not implemented yet, so no "
              f"{REPORT_NAME} was written.")
        print(f"Query the database directly in the meantime, e.g. "
              f"sqlite3 {db_path} '.tables'")

    return 0


def _print_skipped(entries, explicit=False):
    if not entries:
        return
    for entry in entries:
        missing = ", ".join(entry["missing"]) or "no run metadata"
        if explicit:
            print(f"{INFO}INFO:{RESET} Skipped {entry['label'].lower()}: missing {missing}")
        else:
            print(f"  {entry['label']:<24} skipped (missing {missing})")


def run_report_probe(output_dir):
    """Print what a given output directory can support, without building anything."""
    section("DRAKKAR REPORT AVAILABILITY")
    output_path = Path(output_dir)
    if not output_path.is_dir():
        print(f"{ERROR}ERROR:{RESET} Output directory not found: {output_path}")
        return 1
    for entry in probe(output_path):
        status = "available" if entry["available"] else "missing"
        print(f"  {entry['label']:<24} {status}")
        for relative in entry["missing"]:
            print(f"      missing: {relative}")
    return 0
