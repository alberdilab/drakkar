from __future__ import annotations

import re
import sqlite3
import tempfile
import textwrap
import unittest
from datetime import datetime, timedelta, timezone
from pathlib import Path

from drakkar.report import command as report_command
from drakkar.report.render import (
    PAGE_ROWS,
    PALETTE,
    STYLESHEET,
    TREE_MAX_LABELS,
    _runtime_floors,
    render_report,
)
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
    # One row per final bin, with the origin Binette recorded for it: two bins
    # a single binner produced, one both of them produced identically, and one
    # Binette built itself.
    connection.executemany(
        "INSERT INTO assembly_bin (assembly_id, bin_name, is_original, origin) "
        "VALUES (?, ?, ?, ?)",
        [("A1", "A1_bin_1", 1, "metabat"),
         ("A1", "A1_bin_2", 1, "metabat;semibin"),
         ("A1", "A1_bin_3", 0, "binette"),
         ("A1", "A1_bin_4", 1, "semibin"),
         ("A2", "A2_bin_1", 1, "metabat")],
    )
    connection.executemany(
        "INSERT INTO assembly_bin_origin VALUES (?, ?, ?)",
        [("A1", "A1_bin_1", "metabat2"),
         ("A1", "A1_bin_2", "metabat2"),
         ("A1", "A1_bin_2", "semibin2"),
         ("A1", "A1_bin_3", "binette"),
         ("A1", "A1_bin_4", "semibin2"),
         ("A2", "A2_bin_1", "metabat2")],
    )
    log(connection, "assembly", "cataloging", 2)
    log(connection, "assembly_sample", "cataloging", 2)
    log(connection, "assembly_binner", "cataloging", 3)
    log(connection, "assembly_bin", "cataloging", 5)
    log(connection, "assembly_bin_origin", "cataloging", 6)


def seed_dereplication(connection):
    connection.execute(
        "INSERT INTO dereplication VALUES (?, ?, ?, ?, ?, ?, ?)",
        (11, 82.4, 3.1, 0.98, 6, 88.0, 2.4),
    )
    # Four pairs sat close enough for dRep to compare them; the last three bins
    # had no MASH neighbour at all and were never candidates for collapsing.
    clusters = [
        ("G1", "1", "1_1", 1, 0.03), ("G2", "1", "1_1", 0, 0.03),
        ("G3", "2", "2_1", 1, 0.03), ("G4", "2", "2_1", 0, 0.03),
        ("G5", "3", "3_1", 1, 0.05), ("G6", "3", "3_1", 0, 0.05),
        ("G7", "4", "4_1", 1, 0.08), ("G8", "4", "4_2", 1, 0.08),
        ("G9", "5", "5_1", 1, 0.40), ("G10", "6", "6_1", 1, 0.42),
        ("G11", "7", "7_1", 1, 0.44),
    ]
    connection.executemany(
        "INSERT INTO genome_cluster VALUES (?, ?, ?, ?, ?)", clusters
    )
    connection.executemany(
        "INSERT INTO genome_comparison VALUES (?, ?, ?, ?, ?)",
        [("G1", "G2", "1", 0.9990, 0.85),
         ("G3", "G4", "2", 0.9935, 0.81),
         ("G5", "G6", "3", 0.9812, 0.77),
         ("G7", "G8", "4", 0.9702, 0.72)],
    )
    log(connection, "dereplication", "dereplication", 1)
    log(connection, "genome_cluster", "dereplication", len(clusters))
    log(connection, "genome_comparison", "dereplication", 4)


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


def seed_amr(connection):
    # Two assemblies: one where both callers agreed on most loci and geNomad
    # put some of them on a plasmid, one single-caller assembly with a hit the
    # caller reported without coordinates.
    connection.executemany(
        "INSERT INTO amr_assembly VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [("A1", "metagenome", "", 1200, 5_000_000, 5, 4, 2, 6, 2, 3, 0, 1, 3),
         ("A2", "isolate", "Escherichia coli", 80, 4_800_000, 3, 0, 0, 3, 0, 0, 1, 0, 0)],
    )
    loci = [
        ("A1", "L1", "c1", 100, 900, "+", "blaTEM", "blaTEM", "", "", "beta-lactam",
         "", "antibiotic inactivation", "amrfinderplus;rgi", 2, 2,
         "amrfinder_and_rgi", "gene_match"),
        ("A1", "L2", "c1", 2000, 2600, "-", "tetA", "tetA", "", "", "tetracycline",
         "", "antibiotic efflux", "amrfinderplus;rgi", 2, 2,
         "amrfinder_and_rgi", "drug_class_match"),
        ("A1", "L3", "c2", 50, 700, "+", "sul1", "sul1", "", "", "sulfonamide",
         "", "target replacement", "amrfinderplus", 1, 1,
         "amrfinder_only", "single_source"),
        ("A1", "L4", "c2", 3000, 3500, "+", "qnrS", "qnrS", "", "", "quinolone",
         "", "target protection", "rgi", 1, 1, "rgi_only", "single_source"),
        ("A1", "L5", "c3", 10, 400, "+", "blaTEM", "blaTEM", "", "", "beta-lactam",
         "", "antibiotic inactivation", "rgi", 1, 1, "rgi_only", "single_source"),
        ("A1", "L6", "c3", 900, 1400, "+", "mphA", "mphA", "", "", "macrolide",
         "", "antibiotic inactivation", "amrfinderplus", 1, 1,
         "amrfinder_only", "single_source"),
        ("A2", "L7", "d1", 100, 800, "+", "blaTEM", "blaTEM", "", "", "beta-lactam",
         "", "antibiotic inactivation", "amrfinderplus", 1, 1,
         "amrfinder_only", "single_source"),
        ("A2", "L8", "d1", 2000, 2400, "+", "aac(6')-Ib", "aac(6')-Ib", "", "",
         "aminoglycoside", "", "antibiotic inactivation", "amrfinderplus", 1, 1,
         "amrfinder_only", "single_source"),
        ("A2", "L9", "d2", 60, 500, "-", "ermB", "ermB", "", "", "macrolide",
         "", "target alteration", "amrfinderplus", 1, 1,
         "amrfinder_only", "single_source"),
    ]
    connection.executemany(
        "INSERT INTO amr_locus VALUES (" + ", ".join(["?"] * 18) + ")", loci
    )
    hits = [
        ("A1", "H1", "L1", "amrfinderplus", "c1", 100, 900, "+", "blaTEM",
         "beta-lactamase TEM", "", "beta-lactam", "antibiotic inactivation",
         "blaTEM", "EXACTX", "EXACTX", 100.0, 100.0, 540.0, 0),
        ("A1", "H2", "L1", "rgi", "c1", 105, 900, "+", "blaTEM", "TEM-1",
         "ARO:3000873", "beta-lactam", "antibiotic inactivation", "blaTEM",
         "Perfect", "DIAMOND", 99.6, 98.0, 530.0, 0),
        ("A2", "H3", "L7", "amrfinderplus", "d1", 100, 800, "+", "blaTEM",
         "beta-lactamase TEM", "", "beta-lactam", "antibiotic inactivation",
         "blaTEM", "BLASTX", "BLASTX", 97.2, 92.0, 480.0, 1),
    ]
    connection.executemany(
        "INSERT INTO amr_hit VALUES (" + ", ".join(["?"] * 20) + ")", hits
    )
    drugs = [
        ("A1", "H1", "beta-lactam", "L1", "amrfinderplus", "", "antibiotic inactivation", "blaTEM", ""),
        ("A1", "H2", "beta-lactam", "L1", "rgi", "", "antibiotic inactivation", "blaTEM", "ARO:3000873"),
        ("A1", "H4", "tetracycline", "L2", "amrfinderplus", "", "antibiotic efflux", "tetA", ""),
        ("A1", "H5", "sulfonamide", "L3", "amrfinderplus", "", "target replacement", "sul1", ""),
        ("A1", "H6", "quinolone", "L4", "rgi", "", "target protection", "qnr", ""),
        ("A2", "H3", "beta-lactam", "L7", "amrfinderplus", "", "antibiotic inactivation", "blaTEM", ""),
        ("A2", "H7", "macrolide", "L9", "amrfinderplus", "", "target alteration", "ermB", ""),
    ]
    connection.executemany(
        "INSERT INTO amr_drug_class VALUES (" + ", ".join(["?"] * 9) + ")", drugs
    )
    regions = [
        ("A1", "R1", "c1", 1, 40_000, "plasmid", "c1", 40_000, "circular", 0.94,
         0.01, 2.4, 6, 45, "traB;traC", ""),
        ("A1", "R2", "c3", 1, 22_000, "provirus", "c3", 22_000, "linear", 0.81,
         0.04, 1.2, 3, 28, "", "Caudoviricetes"),
    ]
    connection.executemany(
        "INSERT INTO amr_mobility_region VALUES (" + ", ".join(["?"] * 16) + ")",
        regions,
    )
    # L1 sits inside the plasmid, L5 inside the provirus; L2 overlaps both.
    links = [
        ("A1", "L1", "R1", "plasmid", "c1", 800, 1.0, 0.94),
        ("A1", "L2", "R1", "plasmid", "c1", 600, 1.0, 0.94),
        ("A1", "L5", "R2", "provirus", "c3", 390, 1.0, 0.81),
    ]
    connection.executemany(
        "INSERT INTO amr_mobility VALUES (" + ", ".join(["?"] * 8) + ")", links
    )
    log(connection, "amr_assembly", "amr", 2)
    log(connection, "amr_locus", "amr", len(loci))
    log(connection, "amr_hit", "amr", len(hits))
    log(connection, "amr_drug_class", "amr", len(drugs))
    log(connection, "amr_mobility_region", "amr", len(regions))
    log(connection, "amr_mobility", "amr", len(links))


SEEDERS = {
    "preprocessing": seed_preprocessing,
    "cataloging": seed_cataloging,
    "dereplication": seed_dereplication,
    "profiling": seed_profiling,
    "taxonomy": seed_taxonomy,
    "function": seed_function,
    "expression": seed_expression,
    "amr": seed_amr,
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

    def test_ingest_stamps_are_written_for_a_reader_not_a_parser(self):
        _, text, _ = self.render()
        stamped = datetime.now(timezone.utc)
        friendly = f"{stamped.day} {stamped:%b %Y} at {stamped:%H:%M} UTC"
        self.assertIn(friendly, text)
        # The raw ISO stamps the database holds never reach a table cell,
        # neither the ingest log nor the run table. (The embedded Plotly
        # bundle carries ISO strings of its own, hence the cell-only search.)
        cells = re.findall(r"<t[dh][^>]*>([^<]*)</t[dh]>", text)
        self.assertEqual(
            [cell for cell in cells if re.search(r"\d{4}-\d{2}-\d{2}T", cell)], []
        )
        # The run row is seeded at a fixed instant, so its cell is exact.
        self.assertIn("25 Aug 2026 at 10:15 UTC", text)


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

    def test_the_brand_shows_the_logo_above_the_title(self):
        text = self.render()
        brand = text[text.index('<div class="brand">'):text.index("</div>")]
        self.assertIn('<img class="logo" src="data:image/png;base64,', brand)
        self.assertLess(brand.index("<img"), brand.index("<h1>"))
        self.assertIn("<h1>Analysis Report</h1>", brand)

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
        # Preprocessing lists each sample twice: once in the read table and
        # once where the microbial fraction is reported.
        self.assertEqual(text.count("<tr><td>S"), 2 * (PAGE_ROWS + 1))

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

    def test_long_tables_list_every_row_and_are_paged(self):
        """Per-sample listings are not truncated; the browser pages them."""
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        create_schema(connection, drakkar_version="2.1.0")
        seed_preprocessing(connection, samples=150)
        connection.commit()
        connection.close()

        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("preprocessing",))
        text = html_path.read_text(encoding="utf-8")
        self.assertNotIn("Showing the first", text)
        self.assertIn(f'<table data-page-size="{PAGE_ROWS}">', text)
        # Once in the read table, once in the microbial fraction table.
        self.assertEqual(text.count("<tr><td>S"), 300)
        # Sample ids sort as text, so S99 is the last row: it must be present.
        self.assertIn("<td>S99</td>", text)


class PreprocessingLayoutTests(TemporaryRootMixin, unittest.TestCase):
    """Read fates, the gigabases under each mean, and the optional estimates."""

    def section(self, rows):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        create_schema(connection, drakkar_version="2.1.0")
        connection.executemany(
            "INSERT INTO sample (sample_id, reads_pre_fastp, bases_pre_fastp, "
            "reads_post_fastp, host_reads, host_bases, metagenomic_reads, "
            "metagenomic_bases, singlem_fraction, nonpareil_C, "
            "nonpareil_diversity) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        log(connection, "sample", "preprocessing", len(rows))
        connection.commit()
        connection.close()

        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("preprocessing",))
        text = html_path.read_text(encoding="utf-8")
        start = text.index('id="section-preprocessing"')
        return text[start:text.index("</section>", start)]

    def rows(self, estimates=True):
        return [
            (f"S{index}", 1_000_000 * index, 150_000_000 * index,
             950_000 * index, 100_000 * index, 15_000_000 * index,
             850_000 * index, 127_000_000 * index,
             0.8 if estimates else None,
             0.9 if estimates else None,
             18.0 + index if estimates else None)
            for index in (1, 2)
        ]

    def test_mean_reads_carry_the_same_quantity_in_gigabases(self):
        body = self.section(self.rows())
        stats = body[body.index('<div class="stats">'):body.index('<div class="scroll">')]
        # 150 and 300 million bases in, so 0.23 GB on average.
        self.assertIn('<span class="stat-sub">0.23 GB</span>', stats)
        # Host bases run 15 and 30 million, so 0.02 GB on average.
        self.assertIn('<span class="stat-sub">0.02 GB</span>', stats)

    def test_host_reads_are_summarised_before_the_metagenomic_ones(self):
        body = self.section(self.rows())
        self.assertLess(
            body.index("Mean host reads"), body.index("Mean metagenomic reads")
        )

    def test_the_read_fates_chart_is_wide_and_lays_its_legend_in_a_row(self):
        body = self.section(self.rows())
        self.assertIn('<div class="figure wide">', body)
        self.assertIn('"legend":{"title":{"text":""},"orientation":"h"', body)
        # The class only widens the figure if the page carries the rule.
        self.assertIn(".figure.wide", STYLESHEET)

    def test_the_optional_estimates_get_a_subsection_of_their_own(self):
        body = self.section(self.rows())
        self.assertIn("<h3>Microbial fraction and coverage</h3>", body)
        reads = body[:body.index("<h3>Microbial fraction and coverage</h3>")]
        estimates = body[body.index("<h3>Microbial fraction and coverage</h3>"):]
        self.assertNotIn("<th>Microbial fraction</th>", reads)
        # Fractions of one are shown as the percentages they are.
        self.assertIn("<td data-sort=\"0.8\">80.0%</td>", estimates)
        self.assertIn("<th>Nonpareil completeness</th>", estimates)

    def test_the_subsection_is_absent_when_neither_estimate_was_made(self):
        body = self.section(self.rows(estimates=False))
        self.assertIn("<h3>Read fates</h3>", body)
        self.assertNotIn("Microbial fraction and coverage", body)


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

    def test_the_bins_subsection_totals_what_went_in_and_what_came_out(self):
        body = self.section()
        bins = body[body.index("<h3>Bins</h3>"):]
        self.assertIn("Bins produced by the binners", bins)
        # 5 + 6 from A1 and 3 from A2.
        self.assertIn('<span class="stat-value">14</span>', bins)
        self.assertIn("Final bins after Binette", bins)
        # 7 from A1 and 4 from A2.
        self.assertIn('<span class="stat-value">11</span>', bins)
        self.assertIn("78.6% of the bins produced", bins)


class BinFateTests(TemporaryRootMixin, unittest.TestCase):
    """What Binette did with the bins each binner produced."""

    def section(self):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", ("cataloging",))
        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("cataloging",))
        text = html_path.read_text(encoding="utf-8")
        start = text.index('id="section-cataloging"')
        return text[start:text.index("</section>", start)]

    def test_the_fates_follow_the_bins_produced_chart(self):
        body = self.section()
        self.assertLess(
            body.index("<h3>Bins per binner</h3>"),
            body.index("What became of each binner"),
        )

    def test_the_fates_are_drawn_as_a_sankey(self):
        self.assertIn('"type":"sankey"', self.section())

    def test_every_bin_a_binner_produced_is_accounted_for(self):
        body = self.section()
        fates = body[body.index("What became of each binner"):]
        rows = re.findall(r"<tr>(.*?)</tr>", fates, re.S)
        counts = {}
        for row in rows:
            cells = re.findall(r"<td[^>]*>(.*?)</td>", row, re.S)
            if cells:
                counts[cells[0]] = [int(value) for value in cells[1:]]
        # metabat2 produced 8 bins: one kept alone, one shared with semibin2,
        # one more kept alone in A2, and the remaining five replaced.
        self.assertEqual(counts["metabat2"], [8, 3, 1, 5])
        self.assertEqual(counts["semibin2"], [6, 2, 1, 4])
        for produced, kept, _shared, replaced in counts.values():
            self.assertEqual(produced, kept + replaced)

    def test_the_bins_binette_built_itself_are_named(self):
        self.assertIn("Binette built 1 final bin of its own", self.section())

    def test_a_database_without_the_bin_reports_omits_the_subsection(self):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", ("cataloging",))
        connection = sqlite3.connect(db_path)
        connection.execute("DELETE FROM assembly_bin")
        connection.execute("DELETE FROM assembly_bin_origin")
        connection.commit()
        connection.close()
        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("cataloging",))
        text = html_path.read_text(encoding="utf-8")
        self.assertNotIn("What became of each binner", text)
        self.assertIn("<h3>Bins per binner</h3>", text)


class DereplicationLayoutTests(TemporaryRootMixin, unittest.TestCase):
    """Retention, and how the identity threshold acted on the pairs."""

    def section(self, seeded=("dereplication",)):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", seeded)
        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("dereplication",))
        text = html_path.read_text(encoding="utf-8")
        start = text.index('id="section-dereplication"')
        return text[start:text.index("</section>", start)]

    def test_retention_is_highlighted_as_a_share(self):
        body = self.section()
        self.assertIn("Bins retained", body)
        self.assertIn('<span class="stat-value">54.5%</span>', body)
        self.assertIn("5 collapsed into a representative", body)

    def test_the_yield_is_one_stacked_bar_not_two_columns(self):
        body = self.section()
        self.assertIn('"orientation":"h"', body)
        self.assertIn('"barmode":"stack"', body)
        self.assertIn("Collapsed into a representative", body)

    def test_bins_with_no_near_neighbour_are_separated_out(self):
        body = self.section()
        self.assertIn("Bins with a MASH neighbour", body)
        # Eight of the eleven bins sit within 10% MASH distance of another.
        self.assertIn('<span class="stat-value">8</span>', body)
        self.assertIn("72.7% of 11 bins", body)

    def test_pairwise_identities_are_banded_around_the_threshold(self):
        body = self.section()
        bands = body[body.index("Identity band"):]
        self.assertIn("100 – 99.5%", bands)
        self.assertIn("99.5 – 99%", bands)
        # The band between the populated ones is kept, at zero: an empty band
        # beside the threshold is the answer, not a row to hide.
        self.assertIn("99 – 98.5%", bands)
        self.assertIn("97.5 – 97%", bands)

    def test_the_threshold_is_drawn_on_the_histogram(self):
        self.assertIn("98% threshold", self.section())

    def test_a_database_without_the_drep_tables_still_renders_the_summary(self):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", ("dereplication",))
        connection = sqlite3.connect(db_path)
        connection.execute("DELETE FROM genome_cluster")
        connection.execute("DELETE FROM genome_comparison")
        connection.commit()
        connection.close()
        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("dereplication",))
        text = html_path.read_text(encoding="utf-8")
        self.assertIn("Bins retained", text)
        self.assertNotIn("How the identity threshold acted", text)


class TaxonomyLineageTests(TemporaryRootMixin, unittest.TestCase):
    """The per-genome lineage table: one row per bin, one column per rank."""

    LINEAGES = [
        ("MAG_1", "Bacteria", "Bacillota", "Clostridia", "Lachnospirales",
         "Lachnospiraceae", "Blautia", "Blautia wexlerae"),
        ("MAG_2", "Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales",
         "Muribaculaceae", None, None),
        ("MAG_3", "Archaea", "Methanobacteriota", None, None, None, None, None),
    ]

    def section(self):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        create_schema(connection, drakkar_version="2.1.0")
        connection.executemany(
            "INSERT INTO genome_taxonomy (genome_id, domain, phylum, class, "
            '"order", family, genus, species) VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
            self.LINEAGES,
        )
        log(connection, "genome_taxonomy", "taxonomy", len(self.LINEAGES))
        connection.commit()
        connection.close()

        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("taxonomy",))
        text = html_path.read_text(encoding="utf-8")
        start = text.index("<h3>Lineage per genome</h3>")
        return text[start:text.index("</table>", start)]

    def test_every_rank_has_its_own_column(self):
        body = self.section()
        for rank in ("Domain", "Phylum", "Class", "Order", "Family", "Genus",
                     "Species"):
            self.assertIn(f"<th>{rank}</th>", body)
        self.assertIn("<th>Genome</th>", body)

    def test_one_row_per_genome(self):
        body = self.section()
        self.assertEqual(body.count("<tr><td>MAG_"), len(self.LINEAGES))

    def test_ranks_without_a_classification_are_named_not_blanked(self):
        body = self.section()
        # MAG_3 is placed no further than its phylum; the five finer ranks say so.
        row = body[body.index("<tr><td>MAG_3</td>"):]
        row = row[:row.index("</tr>")]
        self.assertEqual(row.count("<td>Unclassified</td>"), 5)

    def test_rows_are_ordered_by_lineage(self):
        # Archaea before Bacteria, so relatives sit together by default.
        body = self.section()
        self.assertLess(body.index("MAG_3"), body.index("MAG_1"))


class GenomeAbundanceTests(TemporaryRootMixin, unittest.TestCase):
    """The per-sample mapping table under the abundance heatmap."""

    def section(self):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", ("preprocessing", "profiling"))
        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("profiling",))
        text = html_path.read_text(encoding="utf-8")
        start = text.index("<h3>Genome abundance</h3>")
        return text[start:text.index("</section>", start)]

    def test_mapped_reads_are_shown_against_what_was_mapped(self):
        body = self.section()
        for header in ("Mapped reads", "Metagenomic reads", "Mapped %"):
            self.assertIn(f"<th>{header}</th>", body)

    def test_the_three_mapping_figures_are_highlighted(self):
        body = self.section()
        for label in ("Mean mapped reads", "Mean metagenomic reads",
                      "Mean mapped %"):
            self.assertIn(f'<span class="stat-label">{label}</span>', body)

    def test_the_metagenomic_columns_drop_out_without_preprocessing(self):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", ("profiling",))
        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("profiling",))
        text = html_path.read_text(encoding="utf-8")
        body = text[text.index("<h3>Genome abundance</h3>"):]
        self.assertIn("<th>Mapped reads</th>", body)
        self.assertNotIn("<th>Metagenomic reads</th>", body)
        self.assertNotIn("<th>Mapped %</th>", body)


class PhylogenyTests(TemporaryRootMixin, unittest.TestCase):
    """The circular phylogeny: the tree, its rings, and what it leaves out."""

    LINEAGES = [
        ("MAG_1", "Bacteria", "Bacillota"),
        ("MAG_2", "Bacteria", "Bacillota"),
        ("MAG_3", "Bacteria", "Bacteroidota"),
        ("MAG_4", "Bacteria", "Bacteroidota"),
    ]
    QUALITY = [
        ("MAG_1", 95.0, 1.0, 3_000_000),
        ("MAG_2", 71.0, 4.0, 1_500_000),
        ("MAG_3", 88.0, 0.5, 2_200_000),
        ("MAG_4", 62.0, 2.5, 4_000_000),
    ]
    TREE = "((MAG_1:0.10,MAG_2:0.12):0.05,(MAG_3:0.08,MAG_4:0.20):0.06);"
    ARCHAEA = "(MAG_5:0.30,MAG_6:0.22);"

    def render(self, tree=TREE, archaea=None, quality=True, lineages=None):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        create_schema(connection, drakkar_version="2.1.0")
        connection.executemany(
            "INSERT INTO genome_taxonomy (genome_id, domain, phylum) "
            "VALUES (?, ?, ?)",
            self.LINEAGES if lineages is None else lineages,
        )
        if quality:
            connection.executemany(
                "INSERT INTO genome (genome_id, completeness, contamination, size) "
                "VALUES (?, ?, ?, ?)",
                self.QUALITY,
            )
            log(connection, "genome", "profiling", len(self.QUALITY))
        if tree is not None:
            connection.execute(
                "INSERT INTO genome_tree VALUES (?, ?, ?)",
                ("bacteria", tree, tree.count(",") + 1),
            )
            log(connection, "genome_tree:bacteria", "taxonomy", 1)
        if archaea is not None:
            connection.execute(
                "INSERT INTO genome_tree VALUES (?, ?, ?)",
                ("archaea", archaea, archaea.count(",") + 1),
            )
            log(connection, "genome_tree:archaea", "taxonomy", 1)
        log(connection, "genome_taxonomy", "taxonomy", len(self.LINEAGES))
        connection.commit()
        connection.close()

        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("taxonomy",))
        return html_path.read_text(encoding="utf-8")

    def figure(self, **kwargs):
        text = self.render(**kwargs)
        start = text.index("<h3>Phylogeny</h3>")
        return text[start:text.index("<h3>", start + 5)]

    def test_the_tree_is_drawn_as_inline_svg(self):
        body = self.figure()
        self.assertIn('<svg class="tree"', body)
        # No second plotting bundle, and nothing to fetch when the page is
        # opened from a file:// URL.
        self.assertNotIn("<img", body)
        self.assertNotIn("<script", body)

    def test_every_tip_carries_the_rings_the_database_can_fill(self):
        body = self.figure()
        for label in ("Phylum", "Genome size", "Completeness", "Contamination"):
            self.assertIn(f">{label}</text>", body)

    def test_a_phylum_keeps_the_colour_the_rest_of_the_section_gives_it(self):
        text = self.render()
        tree = text[text.index("<h3>Phylogeny</h3>"):text.index("<h3>Lineage")]
        counts = text[text.index("<h3>Genomes per phylum</h3>"):]
        # Two genomes each, so the fallback ordering is by genome count and
        # Bacillota leads. Whatever colour it takes, the tree takes too.
        self.assertIn(f'fill="{PALETTE[0]}"', tree)
        bars = counts.index('"marker":{"color":[')
        self.assertIn(PALETTE[0], counts[bars:counts.index("]", bars)])

    def test_each_tip_is_labelled_and_described(self):
        body = self.figure()
        for genome in ("MAG_1", "MAG_2", "MAG_3", "MAG_4"):
            self.assertIn(f">{genome}</text>", body)
        self.assertIn(
            "MAG_1\nPhylum: Bacillota\nGenome size: 3.00 Mbp\n"
            "Completeness: 95.0%\nContamination: 1.00%",
            body,
        )

    def test_the_quality_rings_drop_out_without_the_profiling_tables(self):
        body = self.figure(quality=False)
        self.assertIn(">Phylum</text>", body)
        for label in ("Genome size", "Completeness", "Contamination"):
            self.assertNotIn(f">{label}</text>", body)

    def test_a_genome_the_tree_holds_but_the_tables_do_not_is_still_drawn(self):
        # The tip is placed and named; its rings are blank, and the key says so.
        body = self.figure(
            tree="((MAG_1:0.1,MAG_2:0.1):0.1,MAG_9:0.3);",
            lineages=self.LINEAGES[:2],
        )
        self.assertIn(">MAG_9</text>", body)
        self.assertIn("has no value for it", body)

    def test_a_full_catalogue_is_not_told_to_look_for_blanks(self):
        self.assertNotIn("has no value for it", self.figure())

    def test_bacteria_and_archaea_get_a_figure_each(self):
        body = self.figure(archaea=self.ARCHAEA)
        self.assertEqual(body.count('<svg class="tree"'), 2)
        self.assertIn(">Bacteria</text>", body)
        self.assertIn(">Archaea</text>", body)

    def test_the_key_names_only_the_phyla_the_tree_holds(self):
        body = self.figure()
        key = body[body.index('<div class="tree-key">'):]
        self.assertIn(">Bacillota</span>", key)
        self.assertIn(">Bacteroidota</span>", key)
        self.assertNotIn("Other phyla", key)

    def test_names_are_dropped_once_they_no_longer_fit_round_the_circle(self):
        many = [
            (f"MAG_{index:03d}", "Bacteria", "Bacillota")
            for index in range(1, TREE_MAX_LABELS + 12)
        ]
        tree = "(" + ",".join(f"{row[0]}:0.1" for row in many) + ");"
        body = self.figure(tree=tree, quality=False, lineages=many)
        self.assertIn('<svg class="tree"', body)
        self.assertNotIn(">MAG_001</text>", body)
        # The tips are still there to point at, and still counted.
        self.assertIn("<title>MAG_001", body)
        self.assertIn(f">{len(many)}</text>", body)

    def test_an_unreadable_tree_leaves_the_rest_of_the_section_standing(self):
        text = self.render(tree="((MAG_1:0.1,MAG_2:0.1);")
        self.assertNotIn("<h3>Phylogeny</h3>", text)
        self.assertIn("<h3>Lineage per genome</h3>", text)

    def test_a_taxonomy_without_a_tree_says_nothing_about_phylogeny(self):
        text = self.render(tree=None)
        self.assertNotIn("<h3>Phylogeny</h3>", text)
        self.assertIn("<h3>Genomes per phylum</h3>", text)


class PhylumColourTests(TemporaryRootMixin, unittest.TestCase):
    """Both phylum figures name their colours from one shared assignment."""

    LINEAGES = [
        # Two genomes of the rarer phylum, one of the more abundant one: the
        # two figures therefore rank the phyla in opposite orders.
        ("MAG_1", "Bacteria", "Bacillota"),
        ("MAG_2", "Bacteria", "Bacteroidota"),
        ("MAG_3", "Bacteria", "Bacteroidota"),
    ]
    COUNTS = [("MAG_1", "S1", 9000.0), ("MAG_2", "S1", 500.0),
              ("MAG_3", "S1", 500.0)]

    def section(self):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        create_schema(connection, drakkar_version="2.1.0")
        connection.executemany(
            "INSERT INTO genome_taxonomy (genome_id, domain, phylum) "
            "VALUES (?, ?, ?)", self.LINEAGES,
        )
        connection.executemany(
            "INSERT INTO genome_count (genome_id, sample_id, read_count) "
            "VALUES (?, ?, ?)", self.COUNTS,
        )
        log(connection, "genome_taxonomy", "taxonomy", len(self.LINEAGES))
        connection.commit()
        connection.close()

        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("taxonomy",))
        text = html_path.read_text(encoding="utf-8")
        start = text.index("<h3>Genomes per phylum</h3>")
        return text[start:text.index("</section>", start)]

    def test_a_phylum_keeps_one_colour_across_both_figures(self):
        body = self.section()
        counts, composition = body.split("<h3>Composition per sample</h3>")
        # Bacillota holds the most reads, so it takes the first palette colour
        # in the composition chart; Bacteroidota holds the most genomes, so it
        # leads the counts chart — and must still be drawn in the second.
        self.assertIn(f'{{"marker":{{"color":"{PALETTE[0]}"}},"name":"Bacillota"',
                      composition)
        self.assertIn(
            f'{{"marker":{{"color":"{PALETTE[1]}"}},"name":"Bacteroidota"',
            composition,
        )
        # The counts chart colours its bars by name, in its own bar order.
        self.assertIn(f'"color":["{PALETTE[1]}","{PALETTE[0]}"]', counts)

    def test_the_composition_chart_runs_the_width_of_the_page(self):
        body = self.section()
        composition = body[body.index("<h3>Composition per sample</h3>"):]
        self.assertIn('<div class="figure wide">', composition)
        # A single legend row spread across that width, under the plot.
        self.assertIn('"orientation":"h"', composition)
        self.assertIn('"entrywidthmode":"fraction"', composition)
        # Measured against the whole figure, not the plotting area, so an
        # entry gets as much of the page's width as there is to give it.
        self.assertIn('"xref":"container"', composition)
        self.assertIn('"yref":"container"', composition)


class AmrLayoutTests(TemporaryRootMixin, unittest.TestCase):
    """The AMR section reports loci, not raw caller hits."""

    def section(self):
        root = self.temporary_root()
        db_path = build_db(root / "drakkar.db", ("amr",))
        html_path = root / "drakkar_report.html"
        render_report(db_path, html_path, sections=("amr",))
        text = html_path.read_text(encoding="utf-8")
        start = text.index('id="section-amr"')
        return text[start:text.index("</section>", start)]

    def test_the_section_leads_with_the_reconciled_locus_count(self):
        body = self.section()
        self.assertIn("reconciled into 9 resistance loci", body)
        self.assertIn("Resistance loci", body)

    def test_loci_called_by_both_tools_are_counted_and_shared(self):
        body = self.section()
        self.assertIn("Called by both tools", body)
        # L1 and L2 of the nine seeded loci.
        self.assertIn('<span class="stat-value">2</span>', body)
        self.assertIn("22.2%", body)

    def test_mobile_loci_are_counted_once_per_locus_not_once_per_region(self):
        body = self.section()
        self.assertIn("On a mobile element", body)
        # Three links, three distinct loci; nothing double-counted.
        self.assertIn('<span class="stat-value">3</span>', body)

    def test_caller_support_is_shown_before_the_agreement_between_callers(self):
        body = self.section()
        self.assertIn("<h3>Caller support</h3>", body)
        self.assertLess(
            body.index("<h3>Caller support</h3>"),
            body.index("Agreement on the loci both callers found"),
        )

    def test_the_agreement_labels_are_written_out_not_left_as_slugs(self):
        body = self.section()
        agreement = body[body.index("Agreement on the loci both callers found"):]
        self.assertIn("<td>Same gene</td>", agreement)
        self.assertIn("<td>Same drug class only</td>", agreement)
        self.assertNotIn("gene_match", agreement)

    def test_drug_classes_are_listed_with_the_loci_behind_them(self):
        body = self.section()
        classes = body[body.index("<h3>Drug classes</h3>"):]
        self.assertIn("<td>beta-lactam</td>", classes)
        self.assertIn("<td>macrolide</td>", classes)

    def test_genes_are_ranked_and_carry_their_mobile_count(self):
        body = self.section()
        genes = body[body.index("<h3>Resistance genes</h3>"):]
        self.assertIn("<th>Mobile loci</th>", genes)
        # blaTEM appears in three loci across both assemblies, more than any other.
        first_row = genes[genes.index("<tbody>"):genes.index("</tr>", genes.index("<tbody>"))]
        self.assertIn("blaTEM", first_row)

    def test_the_context_chart_names_the_loci_on_no_mobile_element(self):
        body = self.section()
        context = body[body.index("Mobile-element context"):]
        self.assertIn("Plasmid", context)
        self.assertIn("Provirus", context)
        self.assertIn("No mobile-element call", context)

    def test_hits_without_coordinates_are_reported_rather_than_dropped(self):
        body = self.section()
        self.assertIn("<h3>Hits without coordinates</h3>", body)
        qc = body[body.index("<h3>Hits without coordinates</h3>"):]
        self.assertIn("<td>A1</td>", qc)
        self.assertIn("<td>A2</td>", qc)

    def test_a_run_that_called_nothing_still_reports_what_was_screened(self):
        root = self.temporary_root()
        db_path = root / "drakkar.db"
        connection = connect(db_path)
        create_schema(connection, drakkar_version="2.1.0")
        connection.execute(
            "INSERT INTO amr_assembly VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ("A1", "isolate", "", 40, 4_000_000, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        )
        log(connection, "amr_assembly", "amr", 1)
        connection.commit()
        connection.close()
        html_path = root / "drakkar_report.html"
        outcome = render_report(db_path, html_path, sections=("amr",))
        self.assertEqual(outcome["rendered"], ["amr"])
        text = html_path.read_text(encoding="utf-8")
        self.assertIn("screened 1 assembly", text)
        self.assertIn("0 resistance loci", text)


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

    def test_job_states_are_also_drawn_as_a_stacked_bar(self):
        body = self.section(self.render())
        outcomes = body[body.index("Job outcomes"):body.index("Requested versus used")]
        # One trace per final state, stacked into a single horizontal bar.
        self.assertIn('"barmode":"stack"', outcomes)
        self.assertIn('"orientation":"h"', outcomes)
        for state in ("COMPLETED", "OUT_OF_MEMORY", "No accounting record"):
            self.assertIn(f'"name":"{state}"', outcomes)
        # Success is drawn in the palette's one green, whatever its rank.
        self.assertIn(f'"color":"{PALETTE[2]}"', outcomes)

    def test_rule_highlights_separate_per_job_runtime_from_the_total(self):
        body = self.section(self.render())
        # 1,800 + 5,400 + 3,600 s over the three jobs sacct spoke for, and a
        # fourth launch it had no record of: 10,800 s in total over 4 jobs.
        self.assertIn('<span class="stat-label">Mean runtime per job</span>'
                      '<span class="stat-value">45.0 min</span>', body)
        self.assertIn('<span class="stat-label">Total runtime, all jobs</span>'
                      '<span class="stat-value">3.00 h</span>', body)

    def test_rule_table_and_figure_carry_the_95th_percentile(self):
        body = self.section(self.render())
        self.assertIn("Memory used %, 95th pct", body)
        self.assertIn("Runtime used %, 95th pct", body)
        self.assertIn("Total runtime (h)", body)
        # assembly's two launches used 99.18% and 74.77% of the memory they
        # asked for; the 95th percentile of that pair is 97.96%.
        self.assertIn(">98.0</td>", body)
        # The figure reaches the same ceiling as a whisker off the median bar.
        self.assertIn('"arrayminus":[0.0,0.0]', body)

    def test_the_rule_table_names_the_runtime_floor_it_cannot_go_below(self):
        body = self.section(self.render())
        self.assertIn("Runtime floor (min)", body)
        self.assertIn("Smallest runtime request (min)", body)
        self.assertIn("Jobs at smallest request %", body)
        # assembly's runtime is cap_runtime(max(15, …)) in cataloging.smk, so
        # 15 minutes is the least it can ask for; both its launches asked for
        # 720, well clear of the floor and sized by the input instead.
        self.assertIn('<td data-sort="15">15</td>', body)
        self.assertIn('<td data-sort="720.0">720.0</td>', body)
        self.assertIn('<td data-sort="100.0">100.0</td>', body)

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


class RuntimeFloorTests(unittest.TestCase):
    """The floors the report reads back out of the workflow's rule files."""

    def setUp(self):
        self.floors = _runtime_floors()

    def test_floors_are_read_from_the_rule_definitions(self):
        # cap_runtime(max(15, int(input.size_mb / 20) * …)) in cataloging.smk.
        self.assertEqual(self.floors["assembly"], 15)
        # The annotation rules floor at ten minutes.
        self.assertEqual(self.floors["genomad"], 10)
        self.assertEqual(self.floors["antismash"], 10)

    def test_a_rule_several_modules_define_differently_is_left_out(self):
        # dereplicate is written three times over, at 60, 15 and 15 minutes,
        # and the rule files do not say which module the run used.
        self.assertNotIn("dereplicate", self.floors)

    def test_every_floor_is_a_positive_number_of_minutes(self):
        self.assertTrue(self.floors)
        for rule, floor in self.floors.items():
            with self.subTest(rule=rule):
                self.assertIsInstance(floor, int)
                self.assertGreater(floor, 0)


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

    def only_report(self, root: Path) -> Path:
        """The single rendered report in an output directory."""
        reports = report_command.find_reports(root)
        self.assertEqual(len(reports), 1, f"expected one report, found {reports}")
        return reports[0]

    def test_a_normal_run_writes_both_artifacts(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root), 0)
        self.assertTrue((root / "reporting" / "drakkar.db").exists())
        self.assertTrue(self.only_report(root).exists())

    def test_both_artifacts_land_in_the_reporting_directory(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root), 0)
        self.assertEqual(self.only_report(root).parent, root / "reporting")
        # Nothing is left at the output root, which is the point of the move.
        self.assertFalse((root / "drakkar.db").exists())
        self.assertEqual(list(root.glob(report_command.REPORT_GLOB)), [])

    def test_the_report_name_carries_the_time_it_was_rendered(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root), 0)
        self.assertRegex(
            self.only_report(root).name, r"^drakkar_report_\d{8}-\d{6}\.html$"
        )

    def test_a_later_report_does_not_replace_an_earlier_one(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root), 0)
        first = self.only_report(root)
        # Two runs inside the same second would land on the same name, so the
        # second render is stamped explicitly rather than raced against.
        later = report_command.report_name(
            datetime.now(timezone.utc) + timedelta(seconds=1)
        )
        render_report(
            report_command.database_path(root),
            report_command.reporting_dir(root) / later,
            sections=("preprocessing",),
        )
        self.assertEqual(
            [path.name for path in report_command.find_reports(root)],
            sorted([first.name, later]),
        )

    def test_db_only_still_skips_rendering(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root, db_only=True), 0)
        self.assertEqual(report_command.find_reports(root), [])

    def test_html_only_renders_without_reading_the_source_tables(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root, db_only=True), 0)
        # Removing the source table proves the renderer reads only the database.
        (root / "preprocessing.tsv").unlink()
        self.assertEqual(report_command.run_report(root, html_only=True), 0)
        text = self.only_report(root).read_text(encoding="utf-8")
        self.assertIn('id="section-preprocessing"', text)

    def test_html_only_without_a_database_fails(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root, html_only=True), 1)
        self.assertEqual(report_command.find_reports(root), [])

    def test_db_only_and_html_only_are_mutually_exclusive(self):
        root = self.output_dir()
        self.assertEqual(
            report_command.run_report(root, db_only=True, html_only=True), 1
        )

    def test_html_only_refuses_a_mismatched_schema(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root, db_only=True), 0)
        connection = sqlite3.connect(report_command.database_path(root))
        connection.execute("UPDATE schema_version SET version = 99")
        connection.commit()
        connection.close()
        self.assertEqual(report_command.run_report(root, html_only=True), 1)

    def test_force_rebuilds_the_database_and_renders_again(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root), 0)
        self.assertEqual(report_command.run_report(root, force=True), 0)
        for report in report_command.find_reports(root):
            self.assertIn(
                'id="section-preprocessing"', report.read_text(encoding="utf-8")
            )

    def test_section_selection_reaches_the_page(self):
        root = self.output_dir()
        self.assertEqual(report_command.run_report(root, sections="preprocessing"), 0)
        text = self.only_report(root).read_text(encoding="utf-8")
        self.assertIn('id="section-preprocessing"', text)
        self.assertIn("excluded by", text)


if __name__ == "__main__":
    unittest.main()
