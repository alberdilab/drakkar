from __future__ import annotations

import re
import sqlite3
import tempfile
import textwrap
import unittest
from datetime import datetime, timezone
from pathlib import Path

from drakkar.report import command as report_command
from drakkar.report.render import PAGE_ROWS, TABLE_ROW_LIMIT, render_report
from drakkar.report.schema import SCHEMA_VERSION, connect, create_schema
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
        "assembly_largest_contig, assembly_N50, assembly_gc_percent, "
        "mapping_rate_percent, final_bins, high_quality_bins, medium_quality_bins, "
        "low_quality_bins, bin_mean_completeness, bin_mean_contamination) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [("A1", 1200, 5_000_000, 180_000, 12000, 44.2, 90.0, 7, 3, 3, 1, 88.4, 2.1),
         ("A2", 800, 3_000_000, 95_000, 9000, 41.8, 82.5, 4, 1, 2, 1, 79.6, 3.4)],
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


def seed_resources(connection, benchmark=True):
    connection.execute(
        "INSERT INTO run VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        ("20260825-101500", "2.1.0", "complete", "preprocessing,cataloging",
         "2026-08-25T10:15:00+00:00", "2026-08-25T14:02:11+00:00", "completed",
         "/scratch/run", "drakkar complete"),
    )
    log(connection, "run", "resources", 1)
    if benchmark:
        seed_benchmark(connection)


def seed_benchmark(connection, run_id="20260825-101500", status="generated"):
    """The SLURM accounting projection: roll-up, launches, and rule medians."""
    connection.execute(
        "INSERT INTO run_benchmark VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, "
        "?, ?, ?, ?, ?)",
        (run_id, status, "slurm", None, "2026-08-25T14:03:00+00:00",
         4, 3, 1, 1, 1, 0, 0, 8, 65000.0, 12600.0, 79200.0, 59000.0, 0.7449),
    )
    log(connection, "run_benchmark", "resources", 1)
    if status != "generated":
        return

    jobs = [
        # run, launch, rule, attempt, key, internal, external, wildcards,
        # req cpus, req mem, req runtime, state, exit, alloc cpus, elapsed,
        # cpu time, max rss, cpu eff, mem eff, runtime eff, oom, timeout
        (run_id, 1, "assembly", 1, "assembly|A1", "12", "9001", "assembly=A1",
         8, 65536.0, 720.0, "OUT_OF_MEMORY", "0:125", 8, 1800.0, 12000.0,
         65000.0, 0.833, 0.9918, 0.0417, 1, 0),
        (run_id, 2, "assembly", 2, "assembly|A1", "12", "9002", "assembly=A1",
         8, 131072.0, 720.0, "COMPLETED", "0:0", 8, 5400.0, 38000.0, 98000.0,
         0.8796, 0.7477, 0.125, 0, 0),
        (run_id, 3, "binning", 1, "binning|A1", "14", "9003", "assembly=A1",
         4, 16384.0, 240.0, "COMPLETED", "0:0", 4, 3600.0, 9000.0, 8000.0,
         0.625, 0.4883, 0.25, 0, 0),
        # A launch sacct could not speak for: neither a success nor a failure.
        (run_id, 4, "binning", 1, "binning|A2", "15", "9004", "assembly=A2",
         4, 16384.0, 240.0, None, None, None, None, None, None, None, None,
         None, 0, 0),
    ]
    connection.executemany(
        "INSERT INTO benchmark_job VALUES (" + ", ".join(["?"] * 22) + ")", jobs
    )
    log(connection, "benchmark_job", "resources", len(jobs))

    rules = [
        (run_id, "assembly", 2, 1, 1, 1, 1, 0, 8, 8, 98304.0, 81500.0, 0.8697,
         720.0, 3600.0, 0.0833, 57600.0, 50000.0, 0.868),
        (run_id, "binning", 2, 2, 0, 0, 0, 0, 4, 4, 16384.0, 8000.0, 0.4883,
         240.0, 3600.0, 0.25, 14400.0, 9000.0, 0.625),
    ]
    connection.executemany(
        "INSERT INTO benchmark_rule VALUES (" + ", ".join(["?"] * 19) + ")", rules
    )
    log(connection, "benchmark_rule", "resources", len(rules))


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
        self.assertIn(f"version {SCHEMA_VERSION}", text)
        self.assertIn("20260825-101500", text)
        self.assertIn("Sections rendered", text)

    def test_provenance_lists_the_ingest_log(self):
        _, text, _ = self.render()
        self.assertIn('id="section-provenance"', text)
        self.assertIn("genome_taxonomy", text)


class LayoutTests(TemporaryRootMixin, unittest.TestCase):
    """The page is a sidebar beside one panel at a time, and degrades without JS."""

    def render(self, sections=SECTION_ORDER, seeded=SECTION_ORDER, samples=2):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        try:
            create_schema(connection, drakkar_version="2.1.0")
            for name in seeded:
                if name == "preprocessing":
                    seed_preprocessing(connection, samples=samples)
                else:
                    SEEDERS[name](connection)
            connection.commit()
        finally:
            connection.close()
        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=sections)
        return html_path.read_text(encoding="utf-8")

    def test_the_sidebar_holds_the_details_and_the_contents(self):
        text = self.render()
        sidebar = text[text.index('<aside class="sidebar">'):text.index("</aside>")]
        self.assertIn("Drakkar version", sidebar)
        self.assertIn("Sections rendered", sidebar)
        self.assertIn('<nav class="toc', sidebar)
        self.assertIn('href="#section-taxonomy"', sidebar)
        # Nothing of the sections themselves may leak into it.
        self.assertNotIn("plotly", sidebar)

    def test_each_section_is_a_panel_and_only_the_first_opens(self):
        text = self.render()
        panels = re.findall(r'<section class="panel([^"]*)" id="section-([a-z]+)"', text)
        self.assertEqual([name for _, name in panels],
                         list(SECTION_ORDER) + ["provenance"])
        self.assertEqual([extra for extra, _ in panels if extra], [" is-active"])

    def test_the_first_contents_entry_matches_the_open_panel(self):
        text = self.render(sections=("taxonomy", "resources"))
        self.assertIn('<a href="#section-taxonomy" class="is-active">', text)
        self.assertIn('<section class="panel is-active" id="section-taxonomy"', text)

    def test_the_script_is_inline_and_the_markup_stands_without_it(self):
        text = self.render()
        # No script means no panel switching and no paging: the CSS hides
        # panels only once the script has stamped the document.
        self.assertIn("html.js main .panel { display: none; }", text)
        self.assertIn("document.documentElement.className", text)


class TableTests(TemporaryRootMixin, unittest.TestCase):
    """Averages above every table, and paging past twenty rows."""

    def render(self, samples):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        try:
            create_schema(connection, drakkar_version="2.1.0")
            seed_preprocessing(connection, samples=samples)
            connection.commit()
        finally:
            connection.close()
        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("preprocessing",))
        return html_path.read_text(encoding="utf-8")

    def test_short_tables_are_not_paginated(self):
        text = self.render(samples=PAGE_ROWS)
        self.assertNotIn("<table data-page-size", text)

    def test_long_tables_carry_the_page_size(self):
        text = self.render(samples=PAGE_ROWS + 1)
        self.assertIn(f'<table data-page-size="{PAGE_ROWS}">', text)
        # Paging happens in the browser, so every row is still in the markup.
        self.assertEqual(text.count("<tr><td>S"), PAGE_ROWS + 1)

    def test_column_averages_are_highlighted_above_the_table(self):
        text = self.render(samples=4)
        stats = text[text.index('<div class="stats">'):text.index('<div class="scroll">')]
        self.assertIn("Mean reads in", stats)
        # Reads in run 1000, 2000, 3000, 4000 for the four seeded samples.
        self.assertIn("<span class=\"stat-value\">2,500</span>", stats)
        self.assertIn("Mean metagenomic %", stats)

    def test_a_single_row_is_its_own_average_and_gets_no_highlights(self):
        text = self.render(samples=1)
        self.assertNotIn('<div class="stats">', text)

    def test_numeric_cells_carry_the_value_the_browser_sorts_on(self):
        text = self.render(samples=2)
        # "1,000" would sort as text before "950"; the raw number does not.
        self.assertIn('<td data-sort="1000">1,000</td>', text)
        # Identifiers sort as text and need nothing extra.
        self.assertIn("<td>S1</td>", text)

    def test_headers_and_the_sorting_script_travel_with_the_page(self):
        text = self.render(samples=2)
        # The markup stays plain: the script makes the headers interactive, so
        # no sorting affordance is offered where it would not work.
        self.assertIn("<th>Sample</th>", text)
        self.assertNotIn('<th class="sortable"', text)
        self.assertIn('header.classList.add("sortable")', text)
        self.assertIn('th.sortable[aria-sort="ascending"]::after', text)


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


class CatalogingLayoutTests(TemporaryRootMixin, unittest.TestCase):
    """Assemblies and bins are described by two tables, not one wide one."""

    def section(self):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", ("cataloging",))
        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("cataloging",))
        text = html_path.read_text(encoding="utf-8")
        start = text.index('id="section-cataloging"')
        return text[start:text.index("</section>", start)]

    def test_assembly_and_bin_statistics_are_separate_tables(self):
        body = self.section()
        self.assertIn("<h3>Assemblies</h3>", body)
        self.assertIn("<h3>Bins</h3>", body)
        assemblies = body[body.index("<h3>Assemblies</h3>"):body.index("<h3>Bins</h3>")]
        bins = body[body.index("<h3>Bins</h3>"):]
        self.assertIn("<th>N50</th>", assemblies)
        self.assertIn("<th>Mapping rate %</th>", assemblies)
        self.assertNotIn("<th>Final bins</th>", assemblies)
        self.assertIn("<th>Final bins</th>", bins)
        self.assertIn("<th>Mean completeness</th>", bins)
        self.assertNotIn("<th>N50</th>", bins)

    def test_low_quality_bins_are_not_given_a_column(self):
        # They are discarded downstream, so only the total counts them.
        self.assertNotIn("<th>Low quality</th>", self.section())

    def test_mapping_rates_stay_with_the_assemblies_they_measure(self):
        # They are reads mapped back to the assembly, not to its bins.
        body = self.section()
        self.assertLess(
            body.index("Per-sample mapping rates"), body.index("<h3>Bins</h3>")
        )


class ResourceBenchmarkTests(TemporaryRootMixin, unittest.TestCase):
    """The resources section: job outcomes and requested versus used figures."""

    def render(self, seeder=seed_resources, **kwargs):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        create_schema(connection, drakkar_version="2.1.0")
        seeder(connection, **kwargs)
        connection.commit()
        connection.close()

        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("resources",))
        return html_path.read_text(encoding="utf-8")

    def section(self, text):
        start = text.index('id="section-resources"')
        return text[start:text.index("</section>", start)]

    def test_job_outcomes_are_counted_and_shown_as_percentages(self):
        body = self.section(self.render())
        self.assertIn("Job outcomes", body)
        # Two of four launches completed; one failed; one has no accounting.
        self.assertIn("<span class=\"stat-label\">Successful</span>"
                      "<span class=\"stat-value\">2 (50.0%)</span>", body)
        self.assertIn("<span class=\"stat-label\">Failed</span>"
                      "<span class=\"stat-value\">1 (25.0%)</span>", body)
        self.assertIn("Jobs submitted", body)
        self.assertIn("Relaunches", body)

    def test_final_states_are_broken_out_with_their_share(self):
        body = self.section(self.render())
        self.assertIn("Final state", body)
        self.assertIn("COMPLETED", body)
        self.assertIn("OUT_OF_MEMORY", body)
        self.assertIn("No accounting record", body)
        self.assertIn("ran out of memory", body)

    def test_requested_and_used_resources_appear_per_rule(self):
        body = self.section(self.render())
        self.assertIn("Requested versus used resources, per rule", body)
        self.assertIn("Memory requested (MB)", body)
        self.assertIn("Peak memory (MB)", body)
        self.assertIn("CPU-hours allocated", body)
        # 98,304 MB requested against a 0.8697 median memory efficiency.
        self.assertIn(">98,304</td>", body)
        self.assertIn(">87.0</td>", body)
        # Rules are ranked by the CPU time they consumed.
        self.assertLess(body.index(">assembly<"), body.index(">binning<"))

    def test_requested_and_used_resources_appear_per_job(self):
        body = self.section(self.render())
        self.assertIn("Requested versus used resources, per job", body)
        # The longest launch leads, and carries both sides of the comparison.
        self.assertIn(">131,072</td>", body)
        self.assertIn(">98,000</td>", body)
        self.assertIn("<td>assembly=A1</td>", body)

    def test_seconds_are_reported_as_minutes_against_the_request(self):
        body = self.section(self.render())
        # 5,400 s elapsed against 720 min requested is 90 min, or 12.5%.
        self.assertIn(">90.0</td>", body)
        self.assertIn(">12.5</td>", body)

    def test_a_run_without_usage_figures_explains_itself(self):
        body = self.section(self.render(
            seeder=lambda connection: seed_benchmark(
                connection, status="unsupported_profile"
            ),
        ))
        self.assertIn("slurm profile", body)
        self.assertNotIn("Job outcomes", body)

    def test_the_section_still_renders_without_any_benchmark(self):
        body = self.section(self.render(benchmark=False))
        self.assertIn("20260825-101500", body)
        self.assertNotIn("Job outcomes", body)


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
