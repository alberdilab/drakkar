from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
COUNT_SCRIPT = REPO_ROOT / "drakkar" / "workflow" / "scripts" / "count_total_reads.py"
PROFILING_SCRIPT = REPO_ROOT / "drakkar" / "workflow" / "scripts" / "profiling_genomes_stats.py"
DEREP_SCRIPT = REPO_ROOT / "drakkar" / "workflow" / "scripts" / "dereplicating_stats.py"
MAG_METADATA_SCRIPT = REPO_ROOT / "drakkar" / "workflow" / "scripts" / "mag_metadata.py"
SNAKEFILE = REPO_ROOT / "drakkar" / "workflow" / "Snakefile"


class ProfilingSummaryTests(unittest.TestCase):
    def test_count_total_reads_counts_all_reads_across_fastqs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            r1 = tmpdir_path / "sample_1.fq.gz"
            r2 = tmpdir_path / "sample_2.fq.gz"
            output = tmpdir_path / "sample.totalreads"
            bases_output = tmpdir_path / "sample.totalbases"

            fastq_1 = "@r1\nACGT\n+\n!!!!\n@r2\nTGCA\n+\n!!!!\n"
            fastq_2 = "@r1\nAAAA\n+\n!!!!\n@r2\nCCCC\n+\n!!!!\n"
            with gzip.open(r1, "wt", encoding="utf-8") as handle:
                handle.write(fastq_1)
            with gzip.open(r2, "wt", encoding="utf-8") as handle:
                handle.write(fastq_2)

            subprocess.run(
                [
                    sys.executable,
                    str(COUNT_SCRIPT),
                    str(r1),
                    str(r2),
                    "-o",
                    str(output),
                    "-b",
                    str(bases_output),
                ],
                check=True,
            )

            self.assertEqual(output.read_text(encoding="utf-8").strip(), "4")
            self.assertEqual(bases_output.read_text(encoding="utf-8").strip(), "16")

    def test_profiling_genomes_stats_outputs_mapping_percentage_and_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            mappedreads = tmpdir_path / "sampleA.mappedreads"
            mappedbases = tmpdir_path / "sampleA.mappedbases"
            totalreads = tmpdir_path / "sampleA.totalreads"
            totalbases = tmpdir_path / "sampleA.totalbases"
            output = tmpdir_path / "profiling_genomes.tsv"

            mappedreads.write_text("8\n", encoding="utf-8")
            mappedbases.write_text("1200\n", encoding="utf-8")
            totalreads.write_text("10\n", encoding="utf-8")
            totalbases.write_text("1500\n", encoding="utf-8")

            subprocess.run(
                [
                    sys.executable,
                    str(PROFILING_SCRIPT),
                    "-r", str(mappedreads),
                    "-b", str(mappedbases),
                    "-t", str(totalreads),
                    "-B", str(totalbases),
                    "-o", str(output),
                ],
                check=True,
            )

            df = pd.read_csv(output, sep="\t")
            self.assertEqual(
                list(df.columns),
                ["sample", "input_reads", "input_bases", "reads_mapped", "bases_mapped", "mapping_percentage"],
            )
            self.assertEqual(df.loc[0, "sample"], "sampleA")
            self.assertEqual(df.loc[0, "input_reads"], 10)
            self.assertEqual(df.loc[0, "input_bases"], 1500)
            self.assertEqual(df.loc[0, "reads_mapped"], 8)
            self.assertEqual(float(df.loc[0, "mapping_percentage"]), 80.00)
            self.assertEqual(df.loc[0, "bases_mapped"], 1200)

    def test_dereplicating_stats_summarizes_input_and_output_quality(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            bins_map = tmpdir_path / "bins_to_files.json"
            wdb = tmpdir_path / "Wdb.csv"
            metadata = tmpdir_path / "all_bin_metadata.csv"
            output = tmpdir_path / "dereplicating.tsv"

            bins_map.write_text(
                json.dumps(
                    {
                        "bin1": "/some/path/A.fa.gz",
                        "bin2": "/other/path/B.fna",
                    }
                ),
                encoding="utf-8",
            )
            pd.DataFrame({"genome": ["A.fa"], "score": [100]}).to_csv(wdb, index=False)
            pd.DataFrame(
                {
                    "genome": ["A.fa", "B.fna"],
                    "completeness": [90.0, 80.0],
                    "contamination": [2.0, 4.0],
                }
            ).to_csv(metadata, index=False)

            subprocess.run(
                [
                    sys.executable,
                    str(DEREP_SCRIPT),
                    "--bins-map", str(bins_map),
                    "--wdb", str(wdb),
                    "--ani", "0.95",
                    "--metadata", str(metadata),
                    "-o", str(output),
                ],
                check=True,
            )

            df = pd.read_csv(output, sep="\t")
            self.assertEqual(df.loc[0, "input_bin_number"], 2)
            self.assertEqual(df.loc[0, "input_bin_completeness"], 85.00)
            self.assertEqual(df.loc[0, "input_bin_contamination"], 3.00)
            self.assertEqual(df.loc[0, "dereplication_ani"], 0.95)
            self.assertEqual(df.loc[0, "output_bin_number"], 1)
            self.assertEqual(df.loc[0, "output_bin_completeness"], 90.00)
            self.assertEqual(df.loc[0, "output_bin_contamination"], 2.00)

    def test_mag_metadata_summarizes_dereplicated_genomes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            genome_a = tmpdir_path / "asm1_bin_1.fa"
            genome_b = tmpdir_path / "asm2_bin_3.fa"
            genome_a.write_text(
                ">asm1_bin_1^0\n" + "GC" * 250 + "\n"
                ">asm1_bin_1^1\n" + "AT" * 100 + "\n"
                ">asm1_bin_1^2\n" + "AC" * 50 + "\n",
                encoding="utf-8",
            )
            genome_b.write_text(">asm2_bin_3^0\n" + "ACGT" * 250 + "\n", encoding="utf-8")

            metadata = tmpdir_path / "all_bin_metadata.csv"
            metadata.write_text(
                "genome,completeness,contamination\n"
                "asm1_bin_1.fa,95.5,1.2\n"
                "asm2_bin_3.fa,80.0,4.4\n"
                "asm3_bin_9.fa,70.1,2.0\n",
                encoding="utf-8",
            )

            wdb = tmpdir_path / "Wdb.csv"
            wdb.write_text(
                "genome,cluster,score\n"
                "asm1_bin_1.fa,1_1,92.1\n"
                "asm2_bin_3.fa,2_1,74.5\n",
                encoding="utf-8",
            )

            cdb = tmpdir_path / "Cdb.csv"
            cdb.write_text(
                "genome,secondary_cluster,primary_cluster\n"
                "asm1_bin_1.fa,1_1,1\n"
                "asm3_bin_9.fa,1_1,1\n"
                "asm2_bin_3.fa,2_1,2\n",
                encoding="utf-8",
            )

            output = tmpdir_path / "mags.tsv"
            subprocess.run(
                [
                    sys.executable,
                    str(MAG_METADATA_SCRIPT),
                    "--genomes",
                    str(genome_a),
                    str(genome_b),
                    "--metadata",
                    str(metadata),
                    "--wdb",
                    str(wdb),
                    "--cdb",
                    str(cdb),
                    "-o",
                    str(output),
                ],
                check=True,
            )

            df = pd.read_csv(output, sep="\t").set_index("magid")
            self.assertEqual(list(df.index), ["asm1_bin_1", "asm2_bin_3"])

            first = df.loc["asm1_bin_1"]
            self.assertEqual(first["completeness"], 95.5)
            self.assertEqual(first["contamination"], 1.2)
            self.assertEqual(first["size"], 800)
            self.assertEqual(first["contigs"], 3)
            self.assertEqual(first["longest_contig"], 500)
            self.assertEqual(first["n50"], 500)
            self.assertEqual(first["gc"], 68.75)
            self.assertEqual(first["cluster"], "1_1")
            self.assertEqual(first["cluster_members"], 2)
            self.assertEqual(first["score"], 92.1)

            second = df.loc["asm2_bin_3"]
            self.assertEqual(second["size"], 1000)
            self.assertEqual(second["contigs"], 1)
            self.assertEqual(second["gc"], 50.0)
            self.assertEqual(second["cluster_members"], 1)

    def test_mag_metadata_without_quality_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            genome = tmpdir_path / "asm1_bin_1.fa"
            genome.write_text(">asm1_bin_1^0\n" + "ACGT" * 100 + "\n", encoding="utf-8")

            output = tmpdir_path / "mags.tsv"
            subprocess.run(
                [
                    sys.executable,
                    str(MAG_METADATA_SCRIPT),
                    "--genomes",
                    str(genome),
                    "-o",
                    str(output),
                ],
                check=True,
            )

            df = pd.read_csv(output, sep="\t")
            self.assertEqual(df.loc[0, "magid"], "asm1_bin_1")
            self.assertEqual(df.loc[0, "size"], 400)
            self.assertTrue(pd.isna(df.loc[0, "completeness"]))
            self.assertTrue(pd.isna(df.loc[0, "cluster"]))

    def test_snakefile_targets_root_summary_tables(self) -> None:
        text = SNAKEFILE.read_text(encoding="utf-8")
        self.assertIn('f"{OUTPUT_DIR}/profiling_genomes.tsv"', text)
        self.assertIn('f"{OUTPUT_DIR}/dereplicating.tsv"', text)
        self.assertIn('f"{OUTPUT_DIR}/profiling_genomes/final/mags.tsv"', text)


if __name__ == "__main__":
    unittest.main()
