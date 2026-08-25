from __future__ import annotations

import re
import sqlite3
import tempfile
import textwrap
import unittest
from datetime import datetime, timezone
from pathlib import Path

from drakkar.report import command as report_command
from drakkar.report.render import TABLE_ROW_LIMIT, render_report
from drakkar.report.schema import connect, create_schema
from drakkar.report.sources import SECTION_ORDER


def log(connection, table, section, rows):
    """Stamp ingest_log the way a loader would, so the section counts as present."""
    connection.execute(
        "INSERT OR REPLACE INTO ingest_log "
        "(table_name, section, source_file, source_mtime, source_bytes, "
        " rows_ingested, ingested_at) VALUES (?, ?, ?, ?, ?, ?, ?)",
        (table, section, f"/synthetic/{table}.tsv", 0.0, 0, rows,
         datetime.now(timezone.utc).isoformat()),
    )


def seed_preprocessing(connection, samples=2):
    rows = [
        (f"S{index}", 1000 * index, 950 * index, 100 * index, 850 * index,
         0.8, 18.0 + index)
        for index in range(1, samples + 1)
    ]
    connection.executemany(
        "INSERT INTO sample (sample_id, reads_pre_fastp, reads_post_fastp, "
        "host_reads, metagenomic_reads, singlem_fraction, nonpareil_diversity) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    log(connection, "sample", "preprocessing", len(rows))


def seed_cataloging(connection):
    connection.executemany(
        "INSERT INTO assembly (assembly_id, assembly_contigs, assembly_total_length, "
        "assembly_N50, mapping_rate_percent, final_bins, high_quality_bins) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        [("A1", 1200, 5_000_000, 12000, 90.0, 7, 3),
         ("A2", 800, 3_000_000, 9000, 82.5, 4, 1)],
    )
    connection.executemany(
        "INSERT INTO assembly_sample VALUES (?, ?, ?, ?, ?)",
        [("A1", "S1", 1, 1, 88.5), ("A1", "S2", 1, 1, 91.2)],
    )
    connection.executemany(
        "INSERT INTO assembly_binner VALUES (?, ?, ?)",
        [("A1", "metabat2", 5), ("A1", "semibin2", 6), ("A2", "metabat2", 3)],
    )
    log(connection, "assembly", "cataloging", 2)
    log(connection, "assembly_sample", "cataloging", 2)
    log(connection, "assembly_binner", "cataloging", 3)


def seed_dereplication(connection):
    connection.execute(
        "INSERT INTO dereplication VALUES (?, ?, ?, ?, ?, ?, ?)",
        (11, 82.4, 3.1, 0.98, 6, 88.0, 2.4),
    )
    log(connection, "dereplication", "dereplication", 1)


def seed_profiling(connection, genomes=2):
    genome_rows = [
        (f"MAG_{index}", 90.0 + index, 1.0 + index, 2_400_000, 120, 12000,
         44.2, f"{index}_1", 2, 88.0)
        for index in range(1, genomes + 1)
    ]
    connection.executemany(
        "INSERT INTO genome (genome_id, completeness, contamination, size, contigs, "
        "n50, gc, cluster, cluster_members, score) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        genome_rows,
    )
    counts = [
        (f"MAG_{index}", sample, 1000.0 * index, 20000.0 * index)
        for index in range(1, genomes + 1)
        for sample in ("S1", "S2")
    ]
    connection.executemany("INSERT INTO genome_count VALUES (?, ?, ?, ?)", counts)
    log(connection, "genome", "profiling", len(genome_rows))
    log(connection, "genome_count", "profiling", len(counts))


def seed_taxonomy(connection, genome_ids=("MAG_1", "MAG_2")):
    rows = [
        (genome_id, "Bacteria", "Bacillota", "Clostridia", "Lachnospirales",
         "Lachnospiraceae", "Blautia", "Blautia wexlerae", "topology and ANI",
         "GCF_000159015.1", 97.8, 0.85, None, 92.4, None)
        for genome_id in genome_ids
    ]
    connection.executemany(
        "INSERT INTO genome_taxonomy VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    log(connection, "genome_taxonomy", "taxonomy", len(rows))


def seed_function(connection, genes=3):
    gene_rows = [
        ("MAG_1", f"c1_{index}", "c1", index * 100, index * 100 + 90, "+")
        for index in range(1, genes + 1)
    ]
    connection.executemany("INSERT INTO gene VALUES (?, ?, ?, ?, ?, ?)", gene_rows)
    connection.executemany(
        "INSERT INTO annotation_term VALUES (?, ?, ?, ?)",
        [("kegg", "K00001", "alcohol dehydrogenase", "ko")],
    )
    annotations = [
        ("MAG_1", f"c1_{index}", "kegg", 1, 1, "K00001", 1e-40, 180.2, 1.0)
        for index in range(1, genes)
    ]
    connection.executemany(
        "INSERT INTO gene_annotation VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", annotations
    )
    connection.execute(
        "INSERT INTO cluster_annotation VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        ("MAG_1", "dbcan", "CGC1", "c1", 1, 5000, "CAZyme", 6, "starch", "PUL0012"),
    )
    connection.execute(
        "INSERT INTO annotation_qc VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        ("MAG_1", "gene", "kegg", 500, 420, 80, 0, 410, "merge"),
    )
    log(connection, "gene", "function", len(gene_rows))
    log(connection, "gene_annotation", "function", len(annotations))
    log(connection, "cluster_annotation", "function", 1)
    log(connection, "annotation_qc", "function", 1)


def seed_expression(connection, genes=4, samples=("S1", "S2")):
    gene_rows = [
        (f"MAG_1_c1_{index}", "c1", index * 100, index * 100 + 90, "+", 600)
        for index in range(1, genes + 1)
    ]
    connection.executemany(
        "INSERT INTO expressed_gene VALUES (?, ?, ?, ?, ?, ?)", gene_rows
    )
    counts = [
        (gene[0], sample, float(index * 10))
        for index, gene in enumerate(gene_rows, start=1)
        for sample in samples
    ]
    connection.executemany("INSERT INTO gene_expression VALUES (?, ?, ?)", counts)
    log(connection, "expressed_gene", "expression", len(gene_rows))
    log(connection, "gene_expression", "expression", len(counts))


def seed_resources(connection):
    connection.execute(
        "INSERT INTO run VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        ("20260825-101500", "2.1.0", "complete", "preprocessing,cataloging",
         "2026-08-25T10:15:00+00:00", "2026-08-25T14:02:11+00:00", "completed",
         "/scratch/run", "drakkar complete"),
    )
    log(connection, "run", "resources", 1)


SEEDERS = {
    "preprocessing": seed_preprocessing,
    "cataloging": seed_cataloging,
    "dereplication": seed_dereplication,
    "profiling": seed_profiling,
    "taxonomy": seed_taxonomy,
    "function": seed_function,
    "expression": seed_expression,
    "resources": seed_resources,
}


def build_db(db_path: Path, sections=SECTION_ORDER) -> Path:
    """Create a small synthetic report database holding the named sections."""
    connection = connect(db_path)
    try:
        create_schema(connection, drakkar_version="2.1.0")
        for name in sections:
            SEEDERS[name](connection)
        connection.commit()
    finally:
        connection.close()
    return db_path


class TemporaryRootMixin:
    def temporary_root(self) -> Path:
        """A temporary directory removed when the test ends.

        ``TestCase.enterContext`` would do the same, but it only exists from
        Python 3.11 onwards and Drakkar supports 3.10.
        """
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        return Path(temporary.name)


class RenderTests(TemporaryRootMixin, unittest.TestCase):
    def render(self, sections=SECTION_ORDER, seeded=SECTION_ORDER, root=None):
        root = root or self.temporary_root()
        db_path = build_db(root / "drakkar.db", seeded)
        html_path = root / "drakkar_report.html"
        outcome = render_report(db_path, html_path, sections=sections)
        return html_path, html_path.read_text(encoding="utf-8"), outcome

    def test_report_is_written_next_to_the_database(self):
        html_path, text, _ = self.render()
        self.assertTrue(html_path.exists())
        self.assertTrue(text.startswith("<!DOCTYPE html>"))
        self.assertIn("<title>Drakkar report</title>", text)

    def test_every_seeded_section_is_rendered(self):
        _, text, outcome = self.render()
        self.assertEqual(list(outcome["rendered"]), list(SECTION_ORDER))
        self.assertEqual(outcome["skipped"], [])
        for name in SECTION_ORDER:
            self.assertIn(f'id="section-{name}"', text)

    def test_page_carries_no_external_references(self):
        _, text, _ = self.render()
        # Nothing may be fetched at open time: no external script, stylesheet,
        # image or frame.
        self.assertEqual(re.findall(r"<script[^>]*\ssrc\s*=", text), [])
        self.assertEqual(re.findall(r"<link[^>]*>", text), [])
        self.assertEqual(re.findall(r"<img[^>]+src=[\"']https?:", text), [])
        self.assertNotIn("<iframe", text)
        self.assertIn("<style>", text)

    def test_plotly_bundle_is_embedded_exactly_once(self):
        _, text, _ = self.render()
        self.assertEqual(text.count("plotly.js v"), 1)
        # Several figures share that single bundle.
        self.assertGreater(text.count("plotly-graph-div"), 1)

    def test_summary_header_stamps_versions_and_runs(self):
        _, text, _ = self.render()
        self.assertIn("Report schema", text)
        self.assertIn("version 1", text)
        self.assertIn("20260825-101500", text)
        self.assertIn("Sections rendered", text)

    def test_provenance_lists_the_ingest_log(self):
        _, text, _ = self.render()
        self.assertIn('id="section-provenance"', text)
        self.assertIn("genome_taxonomy", text)


class MissingSectionTests(TemporaryRootMixin, unittest.TestCase):
    def render(self, seeded, sections=SECTION_ORDER):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", seeded)
        html_path = root / "drakkar_report.html"
        outcome = render_report(db_path, html_path, sections=sections)
        return html_path, html_path.read_text(encoding="utf-8"), outcome

    def test_absent_sections_are_omitted_not_raised(self):
        html_path, text, outcome = self.render(seeded=("preprocessing",))
        self.assertTrue(html_path.exists())
        self.assertEqual(outcome["rendered"], ["preprocessing"])
        self.assertIn("taxonomy", outcome["skipped"])
        self.assertNotIn('id="section-taxonomy"', text)

    def test_the_page_names_what_was_unavailable(self):
        _, text, _ = self.render(seeded=("preprocessing",))
        self.assertIn("Sections omitted", text)
        self.assertIn("Functional annotation", text)
        self.assertIn("not present in the report database", text)

    def test_an_empty_database_renders_a_page_rather_than_failing(self):
        _, text, outcome = self.render(seeded=())
        self.assertEqual(outcome["rendered"], [])
        self.assertIn("nothing to show", text)

    def test_deselected_sections_are_reported_separately(self):
        _, text, outcome = self.render(seeded=SECTION_ORDER, sections=("taxonomy",))
        self.assertEqual(outcome["rendered"], ["taxonomy"])
        self.assertEqual(outcome["skipped"], [])
        self.assertIn("preprocessing", outcome["not_selected"])
        self.assertIn("excluded by", text)
        self.assertNotIn('id="section-preprocessing"', text)

    def test_a_section_present_but_empty_is_treated_as_missing(self):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", ("preprocessing",))
        connection = sqlite3.connect(db_path)
        # A loader that found its file but wrote nothing must not produce a
        # section made of empty tables.
        connection.execute("DELETE FROM sample")
        connection.commit()
        connection.close()
        html_path = root / "drakkar_report.html"
        outcome = render_report(db_path, html_path, sections=SECTION_ORDER)
        self.assertEqual(outcome["rendered"], [])
        self.assertIn("preprocessing", outcome["skipped"])

    def test_a_database_missing_a_table_degrades(self):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", ("preprocessing", "taxonomy"))
        connection = sqlite3.connect(db_path)
        connection.execute("DROP TABLE genome_taxonomy")
        connection.commit()
        connection.close()
        outcome = render_report(
            db_path, root / "drakkar_report.html", sections=SECTION_ORDER
        )
        self.assertEqual(outcome["rendered"], ["preprocessing"])
        self.assertIn("taxonomy", outcome["skipped"])


class BoundedOutputTests(TemporaryRootMixin, unittest.TestCase):
    def test_large_tables_are_summarized_not_listed(self):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        create_schema(connection, drakkar_version="2.1.0")
        seed_expression(connection, genes=5000)
        connection.commit()
        connection.close()

        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("expression",))
        text = html_path.read_text(encoding="utf-8")
        # 10,000 expression rows must collapse to one row per sample; no
        # individual gene identifier may reach the page.
        self.assertNotIn("MAG_1_c1_4999", text)
        self.assertIn("Quantification covers 5,000 genes across 2 samples", text)

    def test_long_tables_are_truncated_with_a_note(self):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        create_schema(connection, drakkar_version="2.1.0")
        seed_preprocessing(connection, samples=TABLE_ROW_LIMIT + 50)
        connection.commit()
        connection.close()

        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("preprocessing",))
        text = html_path.read_text(encoding="utf-8")
        self.assertIn(f"Showing the first {TABLE_ROW_LIMIT:,}", text)
        # Sample ids sort as text, so S99 is the last row and must be cut.
        self.assertNotIn("<td>S99</td>", text)


class EscapingTests(TemporaryRootMixin, unittest.TestCase):
    def test_identifiers_are_escaped_into_the_page(self):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        create_schema(connection, drakkar_version="2.1.0")
        seed_taxonomy(connection, genome_ids=("<script>alert(1)</script>",))
        connection.commit()
        connection.close()

        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("taxonomy",))
        text = html_path.read_text(encoding="utf-8")
        self.assertNotIn("<script>alert(1)</script>", text)


class ReportCommandRenderTests(TemporaryRootMixin, unittest.TestCase):
    def output_dir(self) -> Path:
        root = self.temporary_root()
        (root / "preprocessing.tsv").write_text(
            textwrap.dedent("""
                sample\treads_pre_fastp\treads_post_fastp\tmetagenomic_reads
                S1\t1000\t950\t880
                S2\t2000\t1900\t1760
                """).lstrip(),
            encoding="utf-8",
        )
        return root

    def test_a_normal_run_writes_both_artifacts(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root), 0)
        self.assertTrue((root / "drakkar.db").exists())
        self.assertTrue((root / "drakkar_report.html").exists())

    def test_db_only_still_skips_rendering(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root, db_only=True), 0)
        self.assertFalse((root / "drakkar_report.html").exists())

    def test_html_only_renders_without_reading_the_source_tables(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root, db_only=True), 0)
        # Removing the source table proves the renderer reads only the database.
        (root / "preprocessing.tsv").unlink()
        self.assertEqual(report_command.run_report(root, html_only=True), 0)
        text = (root / "drakkar_report.html").read_text(encoding="utf-8")
        self.assertIn('id="section-preprocessing"', text)

    def test_html_only_without_a_database_fails(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root, html_only=True), 1)
        self.assertFalse((root / "drakkar_report.html").exists())

    def test_db_only_and_html_only_are_mutually_exclusive(self):
        root = self.output_dir()
        self.assertEqual(
            report_command.run_report(root, db_only=True, html_only=True), 1
        )

    def test_html_only_refuses_a_mismatched_schema(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root, db_only=True), 0)
        connection = sqlite3.connect(root / "drakkar.db")
        connection.execute("UPDATE schema_version SET version = 99")
        connection.commit()
        connection.close()
        self.assertEqual(report_command.run_report(root, html_only=True), 1)

    def test_force_rebuilds_and_overwrites_the_report(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root), 0)
        first = (root / "drakkar_report.html").read_text(encoding="utf-8")
        self.assertEqual(report_command.run_report(root, force=True), 0)
        second = (root / "drakkar_report.html").read_text(encoding="utf-8")
        self.assertIn('id="section-preprocessing"', first)
        self.assertIn('id="section-preprocessing"', second)

    def test_section_selection_reaches_the_page(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root, sections="preprocessing"), 0)
        text = (root / "drakkar_report.html").read_text(encoding="utf-8")
        self.assertIn('id="section-preprocessing"', text)
        self.assertIn("excluded by", text)


if __name__ == "__main__":
    unittest.main()
