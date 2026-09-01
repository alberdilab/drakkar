"""Report generation: a SQLite projection of a Drakkar output directory, and
the self-contained HTML report rendered from it, both written to the
``reporting/`` subdirectory of the output directory."""

from drakkar.report.command import (
    DATABASE_NAME,
    REPORT_GLOB,
    REPORT_PREFIX,
    REPORT_SUFFIX,
    REPORTING_DIRNAME,
    build_database,
    database_path,
    find_reports,
    report_name,
    reporting_dir,
    run_report,
    run_report_probe,
)
from drakkar.report.schema import SCHEMA_VERSION
from drakkar.report.sources import SECTION_ORDER, parse_sections, probe

__all__ = [
    "DATABASE_NAME",
    "REPORT_GLOB",
    "REPORT_PREFIX",
    "REPORT_SUFFIX",
    "REPORTING_DIRNAME",
    "SCHEMA_VERSION",
    "SECTION_ORDER",
    "build_database",
    "database_path",
    "find_reports",
    "parse_sections",
    "probe",
    "report_name",
    "reporting_dir",
    "run_report",
    "run_report_probe",
]
