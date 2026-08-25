from __future__ import annotations

import lzma
import sqlite3
import tempfile
import textwrap
import unittest
from pathlib import Path

from drakkar.report import command as report_command
from drakkar.report import ingest as report_ingest
from drakkar.report.schema import SCHEMA_VERSION, connect, create_schema
from drakkar.report.sources import SectionError, parse_sections, probe


GENE_COLUMNS = [
    "mag", "gene", "contig", "start", "end", "strand", "source", "method",
    "evidence", "hit_rank", "is_primary", "rank_score", "rank_score_type",
    "annotation_id", "annotation", "annotation_type", "evalue", "bitscore",
    "score", "score_type", "threshold", "identity", "coverage",
    "query_coverage", "target_coverage", "confidence", "alignment_length",
    "query_start", "query_end", "target_start", "target_end", "model_start",
    "model_end", "details",
]


def gene_row(**values):
    row = {column: "" for column in GENE_COLUMNS}
    row["details"] = '{"native": "blob"}'
    row.update({key: str(value) for key, value in values.items()})
    return "\t".join(row[column] for column in GENE_COLUMNS)


def write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(content).lstrip(), encoding="utf-8")


def write_xz(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with lzma.open(path, "wt", encoding="utf-8") as handle:
        handle.write(textwrap.dedent(content).lstrip())


class ReportFixtureMixin:
    def build_output_dir(self, root: Path) -> Path:
        write(root / "preprocessing.tsv", """
            sample\treads_pre_fastp\tbases_pre_fastp\tmetagenomic_reads\tsinglem_fraction\tnonpareil_diversity
            S1\t1000\t150000\t880\t0.85\t18.2
            S2\t2000\t300000\t1960\tNA\t19.1
            """)
        write(root / "cataloging.tsv", """
            assembly\tsamples\tcoverage_samples\tassembly_contigs\tmapping_rate_percent\tsample_mapping_rates\tmetabat2_bins\tmaxbin2_bins\tsemibin2_bins\tcomebin_bins\tfinal_bins
            A1\tS1,S2\tS1,S2\t1200\t90.0\tS1:88.50;S2:91.20\t5\t4\t6\t3\t7
            A2\tS2\tNA\t800\t0\tNA\t3\t2\t4\t2\t4
            """)
        write(root / "dereplicating.tsv", """
            input_bin_number\tinput_bin_completeness\tdereplication_ani\toutput_bin_number
            11\t82.4\t0.98\t6
            """)
        write(root / "profiling_genomes/final/mags.tsv", """
            magid\tcompleteness\tcontamination\tsize\tcontigs\tn50\tgc\tcluster\tcluster_members\tscore
            MAG_A\t95.1\t1.5\t2400000\t120\t12000\t44.2\t1_1\t3\t93.6
            MAG_B\tNA\tNA\t2100000\t180\t9000\tNA\t2_1\t2\tNA
            """)
        write(root / "profiling_genomes/final/counts.tsv", """
            Genome\tS1\tS2
            MAG_A\t15000\t22000
            MAG_B\t8000\t5000
            """)
        write(root / "profiling_genomes/final/bases.tsv", """
            Genome\tS1\tS2
            MAG_A\t2100000\t2350000
            MAG_B\t1700000\t1200000
            """)
        # Raw GTDB-Tk concatenation: packed lineage, blob columns, and the
        # archaeal summary's repeated header row.
        write(root / "annotating/genome_taxonomy.tsv", """
            user_genome\tclassification\tclosest_genome_reference\tclosest_genome_ani\tclosest_genome_af\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\tred_value\twarnings
            MAG_A\td__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia;s__Blautia wexlerae\tGCF_000159015.1\t97.8\t0.85\ttopology and ANI\tN/A\tGCF_001.1, Blautia obeum, 95.0, 88.1, 0.6\t92.4\tNA\t
            user_genome\tclassification\tclosest_genome_reference\tclosest_genome_ani\tclosest_genome_af\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\tred_value\twarnings
            MAG_B\td__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__\tGCF_000016525.1\t98.5\t0.91\ttopology and ANI\tN/A\t\t94.0\t0.87\tGenome has warnings
            """)
        write(root / "annotating/annotation_qc.tsv", """
            mag\tlevel\tsource\treported_records\tretained_records\trejected_records\tunmapped_records\tunique_entities\tfilter_stage
            MAG_A\tgene\tkegg\t500\t420\t80\t0\t410\tmerge
            """)
        write_xz(
            root / "annotating/gene_annotations.tsv.xz",
            "\t".join(GENE_COLUMNS) + "\n" + "\n".join([
                gene_row(mag="MAG_A", gene="c1_1", contig="c1", start=1, end=900,
                         strand="+", source="prodigal", hit_rank=1, is_primary="True",
                         annotation_id="CDS"),
                gene_row(mag="MAG_A", gene="c1_1", contig="c1", start=1, end=900,
                         strand="+", source="kegg", hit_rank=1, is_primary="True",
                         annotation_id="K00001", annotation="alcohol dehydrogenase",
                         annotation_type="ko", evalue="1e-40", bitscore="180.2"),
                gene_row(mag="MAG_A", gene="c1_1", contig="c1", start=1, end=900,
                         strand="+", source="kegg", hit_rank=2, is_primary="False",
                         annotation_id="K00002", annotation="aldehyde reductase",
                         annotation_type="ko", evalue="1e-20", bitscore="95.1"),
                # Negative strand must survive: "-" is a value, not a null marker.
                gene_row(mag="MAG_A", gene="c1_2", contig="c1", start=1000, end=1600,
                         strand="-", source="prodigal", hit_rank=1, is_primary="True",
                         annotation_id="CDS"),
                gene_row(mag="MAG_B", gene="c9_4", contig="c9", start=50, end=800,
                         strand="+", source="kegg", hit_rank=1, is_primary="True",
                         annotation_id="K00001", annotation="alcohol dehydrogenase",
                         annotation_type="ko", evalue="1e-35", bitscore="160.0"),
            ]) + "\n",
        )
        write_xz(root / "annotating/cluster_annotations.tsv.xz", """
            mag\tcluster_id\tcontig\tstart\tend\tsource\tmethod\tevidence\ttype\tgene_count\tsubstrate\tgene_functions\tpul_id\tdetails
            MAG_A\tCGC1\tc1\t1\t5000\tdbcan\trun_dbcan_cgc\tcgc_prediction\tCAZyme\t6\tstarch\tGH13,GT2\tPUL0012\t{}
            """)
        write_xz(root / "expressing/gene_counts.tsv.xz", """
            # Program:featureCounts v2.0.6; Command:featureCounts -F GFF
            Geneid\tChr\tStart\tEnd\tStrand\tLength\tS1\tS2
            MAG_A_c1_1\tc1;c1\t1;400\t300;900\t+;+\t800\t120\t340
            MAG_A_c1_2\tc1\t1000\t1600\t-\t600\t45\t0
            """)
        write(root / "drakkar_20260825-101500.yaml", """
            run_id: '20260825-101500'
            drakkar_version: 2.0.0
            started_at: '2026-08-25T10:15:00+00:00'
            finished_at: '2026-08-25T14:02:11+00:00'
            command: complete
            modules:
            - preprocessing
            - cataloging
            status: completed
            output_directory: /scratch/run
            argv:
            - drakkar
            - complete
            """)
        return root


class SectionParsingTests(unittest.TestCase):
    def test_none_selects_every_section(self):
        self.assertEqual(parse_sections(None)[0], "preprocessing")
        self.assertIn("function", parse_sections(None))

    def test_all_keyword_selects_every_section(self):
        self.assertEqual(parse_sections("all"), parse_sections(None))

    def test_selection_uses_canonical_order(self):
        self.assertEqual(
            parse_sections("function,preprocessing"),
            ("preprocessing", "function"),
        )

    def test_unknown_section_is_rejected_by_name(self):
        with self.assertRaises(SectionError) as raised:
            parse_sections("preprocessing,bogus")
        self.assertIn("bogus", str(raised.exception))

    def test_whitespace_and_case_are_tolerated(self):
        self.assertEqual(parse_sections(" Taxonomy , FUNCTION "), ("taxonomy", "function"))


class ProbeTests(ReportFixtureMixin, unittest.TestCase):
    def test_complete_directory_reports_every_section_available(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            for entry in probe(root):
                self.assertTrue(entry["available"], entry["section"])

    def test_partial_directory_names_the_missing_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root / "preprocessing.tsv", "sample\treads_pre_fastp\nS1\t10\n")
            results = {entry["section"]: entry for entry in probe(root)}
            self.assertTrue(results["preprocessing"]["available"])
            self.assertFalse(results["taxonomy"]["available"])
            self.assertIn(
                "annotating/genome_taxonomy.tsv", results["taxonomy"]["missing"]
            )

    def test_empty_file_does_not_count_as_available(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "preprocessing.tsv").write_text("", encoding="utf-8")
            results = {entry["section"]: entry for entry in probe(root)}
            self.assertFalse(results["preprocessing"]["available"])


class DatabaseBuildTests(ReportFixtureMixin, unittest.TestCase):
    def build(self, root, sections=None, **kwargs):
        db_path = Path(root) / "drakkar.db"
        report_command.build_database(
            root, sections or parse_sections(None), db_path, **kwargs
        )
        connection = sqlite3.connect(db_path)
        connection.row_factory = sqlite3.Row
        return connection

    def test_schema_version_is_stamped_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            rows = connection.execute("SELECT version FROM schema_version").fetchall()
            self.assertEqual([row[0] for row in rows], [SCHEMA_VERSION])
            connection.close()

    def test_gtdbtk_lineage_is_split_into_ranks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT domain, phylum, genus, species, ani, warnings '
                'FROM genome_taxonomy WHERE genome_id = "MAG_A"'
            ).fetchone()
            self.assertEqual(row["domain"], "Bacteria")
            self.assertEqual(row["phylum"], "Bacillota")
            self.assertEqual(row["genus"], "Blautia")
            self.assertEqual(row["species"], "Blautia wexlerae")
            self.assertAlmostEqual(row["ani"], 97.8)
            connection.close()

    def test_repeated_gtdbtk_header_row_is_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            ids = [
                row[0]
                for row in connection.execute("SELECT genome_id FROM genome_taxonomy")
            ]
            self.assertEqual(sorted(ids), ["MAG_A", "MAG_B"])
            connection.close()

    def test_empty_species_rank_becomes_null(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT species FROM genome_taxonomy WHERE genome_id = "MAG_B"'
            ).fetchone()
            self.assertIsNone(row["species"])
            connection.close()

    def test_gene_coordinates_are_stored_once_per_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            genes = connection.execute(
                "SELECT mag, gene, strand FROM gene ORDER BY mag, gene"
            ).fetchall()
            self.assertEqual(
                [(row["mag"], row["gene"]) for row in genes],
                [("MAG_A", "c1_1"), ("MAG_A", "c1_2"), ("MAG_B", "c9_4")],
            )
            connection.close()

    def test_negative_strand_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT strand FROM gene WHERE mag = "MAG_A" AND gene = "c1_2"'
            ).fetchone()
            self.assertEqual(row["strand"], "-")
            connection.close()

    def test_prodigal_rows_are_not_stored_as_annotations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            sources = [
                row[0]
                for row in connection.execute(
                    "SELECT DISTINCT source FROM gene_annotation"
                )
            ]
            self.assertNotIn("prodigal", sources)
            self.assertIn("kegg", sources)
            connection.close()

    def test_annotation_labels_are_deduplicated_into_terms(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            # K00001 is hit by two different genes but stored as one term.
            hits = connection.execute(
                'SELECT COUNT(*) FROM gene_annotation WHERE annotation_id = "K00001"'
            ).fetchone()[0]
            terms = connection.execute(
                'SELECT COUNT(*) FROM annotation_term WHERE annotation_id = "K00001"'
            ).fetchone()[0]
            self.assertEqual(hits, 2)
            self.assertEqual(terms, 1)
            row = connection.execute(
                'SELECT annotation FROM annotation_term WHERE annotation_id = "K00001"'
            ).fetchone()
            self.assertEqual(row["annotation"], "alcohol dehydrogenase")
            connection.close()

    def test_details_column_is_not_carried_into_the_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            columns = {
                row[1]
                for row in connection.execute("PRAGMA table_info(gene_annotation)")
            }
            self.assertNotIn("details", columns)
            connection.close()

    def test_primary_hits_only_drops_secondary_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root, primary_hits_only=True)
            ranks = [
                row[0] for row in connection.execute("SELECT hit_rank FROM gene_annotation")
            ]
            self.assertEqual(set(ranks), {1})
            connection.close()

    def test_secondary_hits_are_kept_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            ranks = sorted(
                row[0] for row in connection.execute("SELECT hit_rank FROM gene_annotation")
            )
            self.assertIn(2, ranks)
            connection.close()

    def test_packed_sample_mapping_rates_are_normalized(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            rows = connection.execute(
                'SELECT sample_id, mapping_rate_percent FROM assembly_sample '
                'WHERE assembly_id = "A1" ORDER BY sample_id'
            ).fetchall()
            self.assertEqual(
                [(row["sample_id"], row["mapping_rate_percent"]) for row in rows],
                [("S1", 88.5), ("S2", 91.2)],
            )
            connection.close()

    def test_binner_counts_are_normalized(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            rows = dict(
                connection.execute(
                    'SELECT binner, bin_count FROM assembly_binner WHERE assembly_id = "A1"'
                ).fetchall()
            )
            self.assertEqual(
                rows, {"metabat2": 5, "maxbin2": 4, "semibin2": 6, "comebin": 3}
            )
            connection.close()

    def test_counts_and_bases_are_melted_into_long_form(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT read_count, covered_bases FROM genome_count '
                'WHERE genome_id = "MAG_A" AND sample_id = "S2"'
            ).fetchone()
            self.assertEqual(row["read_count"], 22000.0)
            self.assertEqual(row["covered_bases"], 2350000.0)
            connection.close()

    def test_na_values_become_null_not_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT completeness, gc FROM genome WHERE genome_id = "MAG_B"'
            ).fetchone()
            self.assertIsNone(row["completeness"])
            self.assertIsNone(row["gc"])
            connection.close()

    def test_featurecounts_comment_line_and_packed_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT contig, start, end, strand, length FROM expressed_gene '
                'WHERE gene_id = "MAG_A_c1_1"'
            ).fetchone()
            # "1;400" / "300;900" collapse to the meta-feature bounds.
            self.assertEqual((row["contig"], row["start"], row["end"]), ("c1", 1, 900))
            self.assertEqual(row["strand"], "+")
            self.assertEqual(row["length"], 800)
            counts = dict(
                connection.execute(
                    'SELECT sample_id, count FROM gene_expression '
                    'WHERE gene_id = "MAG_A_c1_1"'
                ).fetchall()
            )
            self.assertEqual(counts, {"S1": 120.0, "S2": 340.0})
            connection.close()

    def test_run_metadata_is_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute("SELECT * FROM run").fetchone()
            self.assertEqual(row["run_id"], "20260825-101500")
            self.assertEqual(row["command"], "complete")
            self.assertEqual(row["modules"], "preprocessing,cataloging")
            self.assertEqual(row["status"], "completed")
            connection.close()

    def test_ingest_log_records_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT section, source_file, rows_ingested FROM ingest_log '
                'WHERE table_name = "genome_taxonomy"'
            ).fetchone()
            self.assertEqual(row["section"], "taxonomy")
            self.assertTrue(row["source_file"].endswith("genome_taxonomy.tsv"))
            self.assertEqual(row["rows_ingested"], 2)
            connection.close()

    def test_selected_sections_leave_other_tables_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root, sections=parse_sections("taxonomy"))
            taxonomy = connection.execute(
                "SELECT COUNT(*) FROM genome_taxonomy"
            ).fetchone()[0]
            samples = connection.execute("SELECT COUNT(*) FROM sample").fetchone()[0]
            self.assertEqual(taxonomy, 2)
            self.assertEqual(samples, 0)
            connection.close()

    def test_chunked_ingest_matches_single_pass_ingest(self):
        original = report_ingest.CHUNK_SIZE
        report_ingest.CHUNK_SIZE = 2
        try:
            with tempfile.TemporaryDirectory() as tmp:
                root = self.build_output_dir(Path(tmp))
                connection = self.build(root)
                annotations = connection.execute(
                    "SELECT COUNT(*) FROM gene_annotation"
                ).fetchone()[0]
                genes = connection.execute("SELECT COUNT(*) FROM gene").fetchone()[0]
                terms = connection.execute(
                    "SELECT COUNT(*) FROM annotation_term"
                ).fetchone()[0]
                self.assertEqual(annotations, 3)
                self.assertEqual(genes, 3)
                self.assertEqual(terms, 2)
                connection.close()
        finally:
            report_ingest.CHUNK_SIZE = original


class MissingSourceTests(unittest.TestCase):
    def test_loaders_return_none_when_sources_are_absent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            connection = connect(root / "drakkar.db")
            create_schema(connection)
            for name, loader in report_ingest.SECTION_LOADERS.items():
                self.assertIsNone(loader(connection, root), name)
            connection.close()


class ReportCommandTests(ReportFixtureMixin, unittest.TestCase):
    def test_missing_output_directory_fails(self):
        self.assertEqual(report_command.run_report("/nonexistent-drakkar-dir"), 1)

    def test_unknown_section_fails_without_building(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            self.assertEqual(report_command.run_report(root, sections="bogus"), 1)
            self.assertFalse((root / "drakkar.db").exists())

    def test_empty_directory_fails_with_no_sections(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(report_command.run_report(tmp, db_only=True), 1)

    def test_successful_build_writes_the_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            self.assertEqual(report_command.run_report(root, db_only=True), 0)
            self.assertTrue((root / "drakkar.db").exists())

    def test_schema_mismatch_blocks_until_forced(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            self.assertEqual(report_command.run_report(root, db_only=True), 0)
            connection = sqlite3.connect(root / "drakkar.db")
            connection.execute("UPDATE schema_version SET version = 99")
            connection.commit()
            connection.close()
            self.assertEqual(report_command.run_report(root, db_only=True), 1)
            self.assertEqual(
                report_command.run_report(root, db_only=True, force=True), 0
            )

    def test_rebuild_is_idempotent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            report_command.run_report(root, db_only=True)
            report_command.run_report(root, db_only=True, force=True)
            connection = sqlite3.connect(root / "drakkar.db")
            count = connection.execute("SELECT COUNT(*) FROM gene_annotation").fetchone()[0]
            self.assertEqual(count, 3)
            connection.close()


if __name__ == "__main__":
    unittest.main()
