from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "drakkar" / "workflow" / "Snakefile"
CATALOGING_RULES = ROOT / "drakkar" / "workflow" / "rules" / "cataloging.smk"
CONTIG_TO_BIN_SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "contig_to_bin.py"


def run_script(*args: str) -> None:
    subprocess.run([sys.executable, str(CONTIG_TO_BIN_SCRIPT), *args], check=True)


class ContigToBinTests(unittest.TestCase):
    def test_cataloging_workflow_exports_the_contig_to_bin_tables(self) -> None:
        snakefile = SNAKEFILE.read_text(encoding="utf-8")
        rules = CATALOGING_RULES.read_text(encoding="utf-8")

        self.assertIn('f"{OUTPUT_DIR}/cataloging/final/all_contig_to_bin.csv"', snakefile)
        self.assertIn("rule contig_to_bin:", rules)
        self.assertIn("rule all_contig_to_bin:", rules)
        self.assertIn("contig_to_bin.py", rules)
        self.assertIn("bins=get_final_bin_fastas", rules)

    def test_table_maps_bin_contigs_to_their_assembly(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            bin1 = tmp_path / "assembly1_bin_1.fa"
            bin2 = tmp_path / "assembly1_bin_2.fa"
            # Megahit contigs are renamed <assembly>_<n> by the assembly rule and
            # binette keeps those names, so they are the link to the assembly.
            bin1.write_text(
                ">assembly1_1 flag=1 multi=3.0\nACGT\n>assembly1_7\nTTTT\n",
                encoding="utf-8",
            )
            bin2.write_text(">assembly1_3\nGGGG\n", encoding="utf-8")
            table = tmp_path / "assembly1.tsv"

            run_script("--assembly", "assembly1", "-o", str(table), str(bin1), str(bin2))

            rows = list(csv.DictReader(table.open(encoding="utf-8"), delimiter="\t"))
            self.assertEqual(
                [(row["contig"], row["bin"], row["assembly"]) for row in rows],
                [
                    ("assembly1_1", "assembly1_bin_1", "assembly1"),
                    ("assembly1_7", "assembly1_bin_1", "assembly1"),
                    ("assembly1_3", "assembly1_bin_2", "assembly1"),
                ],
            )

    def test_assemblies_without_bins_export_a_header_only_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            table = Path(tmpdir) / "assembly2.tsv"

            run_script("--assembly", "assembly2", "-o", str(table))

            self.assertEqual(
                table.read_text(encoding="utf-8"), "contig\tbin\tassembly\n"
            )

    def test_merge_combines_per_assembly_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            first = tmp_path / "assembly1.tsv"
            first.write_text(
                "contig\tbin\tassembly\nassembly1_1\tassembly1_bin_1\tassembly1\n",
                encoding="utf-8",
            )
            second = tmp_path / "assembly2.tsv"
            second.write_text("contig\tbin\tassembly\n", encoding="utf-8")
            third = tmp_path / "assembly3.tsv"
            third.write_text(
                "contig\tbin\tassembly\nassembly3_5\tassembly3_bin_2\tassembly3\n",
                encoding="utf-8",
            )
            merged = tmp_path / "all_contig_to_bin.csv"

            run_script("--merge", "-o", str(merged), str(first), str(second), str(third))

            rows = list(csv.DictReader(merged.open(encoding="utf-8")))
            self.assertEqual(
                [(row["contig"], row["bin"], row["assembly"]) for row in rows],
                [
                    ("assembly1_1", "assembly1_bin_1", "assembly1"),
                    ("assembly3_5", "assembly3_bin_2", "assembly3"),
                ],
            )


if __name__ == "__main__":
    unittest.main()
