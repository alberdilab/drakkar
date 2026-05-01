from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "drakkar" / "workflow" / "config.yaml"
SNAKEFILE = ROOT / "drakkar" / "workflow" / "Snakefile"
CATALOGING_RULES = ROOT / "drakkar" / "workflow" / "rules" / "cataloging.smk"
CATALOGING_STATS_SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "cataloging_stats.py"


class CatalogingStatsTests(unittest.TestCase):
    def test_cataloging_workflow_targets_quast_flagstat_and_cataloging_table(self) -> None:
        config = CONFIG_PATH.read_text(encoding="utf-8")
        snakefile = SNAKEFILE.read_text(encoding="utf-8")
        rules = CATALOGING_RULES.read_text(encoding="utf-8")

        self.assertIn('QUAST_MODULE: "quast/5.3.0"', config)
        self.assertIn('f"{OUTPUT_DIR}/cataloging.tsv"', snakefile)
        self.assertIn("rule assembly_quast:", rules)
        self.assertIn("quast.py", rules)
        self.assertIn("rule assembly_flagstat:", rules)
        self.assertIn("samtools flagstat", rules)
        self.assertIn("rule cataloging_stats:", rules)
        self.assertIn("cataloging_stats.py", rules)

    def test_cataloging_stats_script_writes_assembly_mapping_and_binning_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            data_dir = tmp_path / "data"
            quast_dir = tmp_path / "cataloging" / "quast" / "assembly1"
            bowtie_dir = tmp_path / "cataloging" / "bowtie2" / "assembly1"
            final_dir = tmp_path / "cataloging" / "final"
            for directory in [data_dir, quast_dir, bowtie_dir, final_dir]:
                directory.mkdir(parents=True)

            assembly_to_samples = data_dir / "assembly_to_samples.json"
            assembly_to_samples.write_text(
                json.dumps({"assembly1": ["sample1", "sample2"]}),
                encoding="utf-8",
            )

            quast_report = quast_dir / "report.tsv"
            quast_report.write_text(
                "Assembly\tassembly1\n"
                "# contigs\t12\n"
                "Largest contig\t5000\n"
                "Total length\t42000\n"
                "GC (%)\t48.7\n"
                "N50\t3500\n"
                "N75\t2200\n"
                "L50\t4\n"
                "L75\t8\n",
                encoding="utf-8",
            )

            flagstat1 = bowtie_dir / "sample1.flagstat.txt"
            flagstat1.write_text(
                "100 + 0 in total (QC-passed reads + QC-failed reads)\n"
                "80 + 0 mapped (80.00% : N/A)\n",
                encoding="utf-8",
            )
            flagstat2 = bowtie_dir / "sample2.flagstat.txt"
            flagstat2.write_text(
                "50 + 0 in total (QC-passed reads + QC-failed reads)\n"
                "25 + 0 mapped (50.00% : N/A)\n",
                encoding="utf-8",
            )

            metabat2 = tmp_path / "cataloging" / "metabat2" / "assembly1" / "assembly1.tsv"
            maxbin2 = tmp_path / "cataloging" / "maxbin2" / "assembly1" / "assembly1.tsv"
            semibin2 = tmp_path / "cataloging" / "semibin2" / "assembly1" / "assembly1.tsv"
            for path, rows in [(metabat2, 2), (maxbin2, 1), (semibin2, 3)]:
                path.parent.mkdir(parents=True)
                path.write_text("\n".join(f"c{i}\tb{i}" for i in range(rows)) + "\n", encoding="utf-8")

            bins = final_dir / "assembly1.tsv"
            bins.write_text(
                "bin_id\tcompleteness\tcontamination\tscore\tsize\tN50\tcontig_count\n"
                "1\t95\t2\t90\t10000\t4000\t10\n"
                "2\t60\t8\t50\t5000\t1500\t6\n"
                "3\t30\t12\t10\t2000\t800\t4\n",
                encoding="utf-8",
            )

            output = tmp_path / "cataloging.tsv"
            subprocess.run(
                [
                    sys.executable,
                    str(CATALOGING_STATS_SCRIPT),
                    "--assembly-to-samples",
                    str(assembly_to_samples),
                    "--quast",
                    str(quast_report),
                    "--flagstat",
                    str(flagstat1),
                    str(flagstat2),
                    "--metabat2",
                    str(metabat2),
                    "--maxbin2",
                    str(maxbin2),
                    "--semibin2",
                    str(semibin2),
                    "--bins",
                    str(bins),
                    "-o",
                    str(output),
                ],
                check=True,
            )

            table = pd.read_csv(output, sep="\t")
            row = table.iloc[0].to_dict()
            self.assertEqual(row["assembly"], "assembly1")
            self.assertEqual(row["samples"], "sample1,sample2")
            self.assertEqual(row["coverage_samples"], "sample1,sample2")
            self.assertEqual(row["assembly_contigs"], 12)
            self.assertEqual(row["assembly_total_length"], 42000)
            self.assertEqual(row["mapped_reads"], 105)
            self.assertEqual(row["total_reads"], 150)
            self.assertEqual(row["mapping_rate_percent"], 70)
            self.assertEqual(row["sample_mapping_rates"], "sample1:80.00;sample2:50.00")
            self.assertEqual(row["metabat2_bins"], 2)
            self.assertEqual(row["maxbin2_bins"], 1)
            self.assertEqual(row["semibin2_bins"], 3)
            self.assertEqual(row["final_bins"], 3)
            self.assertEqual(row["high_quality_bins"], 1)
            self.assertEqual(row["medium_quality_bins"], 1)
            self.assertEqual(row["low_quality_bins"], 1)
            self.assertEqual(row["best_bin"], 1)
            self.assertEqual(row["best_bin_score"], 90)


if __name__ == "__main__":
    unittest.main()
