"""The `drakkar reporting` command.

Builds ``drakkar.db``, a lean SQLite projection of whatever a Drakkar output
directory happens to contain, and renders ``drakkar_report_<timestamp>.html``
from it. The rendering step reads the database and nothing else, so the
database is always the single source of truth and the source tables are read
exactly once.

Both artefacts live in a ``reporting/`` subdirectory of the output directory,
alongside the ``preprocessing/``, ``cataloging/`` and ``annotating/``
directories the workflows write, so the output root stays a short list of
summary tables rather than accumulating one file per render.
"""

from datetime import datetime, timezone
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

# The reporting artefacts share one directory of their own, named after the
# command that writes them and after the gerund convention the workflow output
# directories already follow.
REPORTING_DIRNAME = "reporting"

DATABASE_NAME = "drakkar.db"

# The rendered report is stamped with the time it was rendered rather than
# given one fixed name, so that re-reporting an output directory — after a
# further workflow has run, or with a different --sections selection — keeps
# the earlier report instead of overwriting it. The stamp is UTC and uses the
# same format as a run id, so a report sorts next to the run it describes.
REPORT_PREFIX = "drakkar_report_"
REPORT_SUFFIX = ".html"
REPORT_GLOB = f"{REPORT_PREFIX}*{REPORT_SUFFIX}"


def report_name(timestamp=None):
    """Return the file name for a report rendered at ``timestamp`` (now by default)."""
    stamp = timestamp if timestamp is not None else datetime.now(timezone.utc)
    return f"{REPORT_PREFIX}{stamp.strftime('%Y%m%d-%H%M%S')}{REPORT_SUFFIX}"


def reporting_dir(output_dir):
    """Return the directory holding the report database and the rendered pages."""
    return Path(output_dir) / REPORTING_DIRNAME


def database_path(output_dir):
    """Return the path the report database is written to."""
    return reporting_dir(output_dir) / DATABASE_NAME


def legacy_database_path(output_dir):
    """Return where builds before the reporting directory put the database."""
    return Path(output_dir) / DATABASE_NAME


def find_reports(output_dir):
    """Return the rendered reports in an output directory, oldest first.

    Reports rendered before the reporting directory existed sat at the output
    root, so both places are searched and the timestamped names — not the
    directories — decide the order.
    """
    output_path = Path(output_dir)
    found = []
    for directory in (reporting_dir(output_path), output_path):
        try:
            found.extend(directory.glob(REPORT_GLOB))
        except OSError:
            continue
    return sorted(found, key=lambda path: path.name)


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


def _adopt_legacy_database(output_path, db_path):
    """Move a pre-``reporting/`` database into the reporting directory.

    Databases built before the reporting directory existed sit at the output
    root. Moving one — with its WAL sidecars, which hold committed pages that
    have not been checkpointed yet — keeps it reusable instead of rebuilding a
    second copy beside it. If the move cannot be made, the database is read and
    written where it already is rather than the run failing over tidiness.
    """
    legacy = legacy_database_path(output_path)
    if db_path.exists() or not legacy.exists():
        return db_path
    try:
        db_path.parent.mkdir(parents=True, exist_ok=True)
        for suffix in ("-wal", "-shm"):
            sidecar = Path(str(legacy) + suffix)
            if sidecar.exists():
                sidecar.replace(Path(str(db_path) + suffix))
        legacy.replace(db_path)
    except OSError:
        return legacy
    print(f"{INFO}INFO:{RESET} Moved the existing report database into "
          f"{REPORTING_DIRNAME}/: {db_path}")
    return db_path


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
        print("Rebuild it from the source tables with 'drakkar reporting --force'.")
        return False
    return True


def run_report(output_dir, sections=None, db_only=False, html_only=False,
               force=False, primary_hits_only=False):
    """Entry point for `drakkar reporting`."""
    section("BUILDING DRAKKAR REPORT")

    output_path = Path(output_dir)
    if not output_path.is_dir():
        print(f"{ERROR}ERROR:{RESET} Output directory not found: {output_path}")
        return 1

    if db_only and html_only:
        print(f"{ERROR}ERROR:{RESET} --db-only and --html-only are mutually exclusive.")
        return 1

    try:
        requested = parse_sections(sections)
    except SectionError as exc:
        print(f"{ERROR}ERROR:{RESET} {exc}")
        return 1

    db_path = _adopt_legacy_database(output_path, database_path(output_path))

    if html_only:
        if not _check_renderable_database(db_path):
            return 1
        print(f"{INFO}INFO:{RESET} Re-rendering from the existing database: {db_path}")
    else:
        if not _check_existing_database(db_path, force):
            return 1

        try:
            db_path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            print(f"{ERROR}ERROR:{RESET} Cannot create the reporting directory: "
                  f"{db_path.parent}")
            print(f"{exc.__class__.__name__}: {exc}")
            return 1

        results = build_database(
            output_path, requested, db_path, primary_hits_only=primary_hits_only
        )

        ingested = [entry for entry in results if entry["available"]]
        missing = [entry for entry in results if not entry["available"]]

        if not ingested:
            print(f"{ERROR}ERROR:{RESET} No report sections could be built from: {output_path}")
            print("None of the expected Drakkar output tables were found there.")
            _print_skipped(missing)
            return 1

        print(f"{INFO}INFO:{RESET} Report database written: {db_path}")
        for entry in ingested:
            rows = entry["rows"]
            detail = f"{rows} rows" if rows is not None else "no rows"
            print(f"  {entry['label']:<24} {detail}")

        # Sections the user explicitly asked for but whose inputs are absent are
        # named individually; an unfiltered run just lists what was skipped.
        _print_skipped(missing, explicit=sections is not None)

        if db_only:
            return 0

    return _render_html(db_path, reporting_dir(output_path) / report_name(), requested)


def _check_renderable_database(db_path):
    """Verify that --html-only has a database of the expected schema to read."""
    if not db_path.exists():
        print(f"{ERROR}ERROR:{RESET} No report database to render: {db_path}")
        print("Build it first with 'drakkar reporting' (without --html-only).")
        return False
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
        print("Rebuild it from the source tables with 'drakkar reporting --force'.")
        return False
    return True


def _render_html(db_path, html_path, requested):
    """Render the HTML report, reporting what the database could not supply."""
    # Imported here so that the plotly import cost is paid only by a render.
    from drakkar.report.render import render_report

    try:
        outcome = render_report(db_path, html_path, sections=requested)
    except Exception as exc:
        print(f"{ERROR}ERROR:{RESET} Could not render {html_path.name}: "
              f"{exc.__class__.__name__}: {exc}")
        return 1

    if not outcome["rendered"]:
        print(f"{ERROR}ERROR:{RESET} The report database holds none of the requested "
              f"sections, so {html_path.name} would be empty.")
        return 1

    print(f"{INFO}INFO:{RESET} HTML report written: {html_path}")
    for name in outcome["rendered"]:
        print(f"  {SECTION_SOURCES[name]['label']:<24} rendered")
    for name in outcome["skipped"]:
        label = SECTION_SOURCES[name]["label"]
        print(f"  {label:<24} not rendered (absent from the database)")
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
