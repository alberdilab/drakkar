"""Report generation: a lean SQLite projection of a Drakkar output directory."""

from drakkar.report.command import (
    DATABASE_NAME,
    REPORT_NAME,
    build_database,
    run_report,
    run_report_probe,
)
from drakkar.report.schema import SCHEMA_VERSION
from drakkar.report.sources import SECTION_ORDER, parse_sections, probe

__all__ = [
    "DATABASE_NAME",
    "REPORT_NAME",
    "SCHEMA_VERSION",
    "SECTION_ORDER",
    "build_database",
    "parse_sections",
    "probe",
    "run_report",
    "run_report_probe",
]
