from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from drakkar import cli as cli_module


class QualityFileTests(unittest.TestCase):
    def test_validate_and_write_quality_file_accepts_extensionless_genome_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            data_dir = tmp_path / "data"
            data_dir.mkdir(parents=True)
            (data_dir / "bins_to_files.json").write_text(
                json.dumps({"MAG_01": "/tmp/MAG_01.fna.gz"}),
                encoding="utf-8",
            )

            quality_file = tmp_path / "mag_qualities.csv"
            quality_file.write_text(
                "genome,completeness,contamination\n"
                "MAG_01,97.3,2.1\n",
                encoding="utf-8",
            )

            self.assertTrue(cli_module.validate_and_write_quality_file(str(quality_file), str(tmp_path)))

            output = pd.read_csv(tmp_path / "cataloging" / "final" / "all_bin_metadata.csv")
            self.assertEqual(output.loc[0, "genome"], "MAG_01.fna")
            self.assertAlmostEqual(output.loc[0, "completeness"], 97.3)
            self.assertAlmostEqual(output.loc[0, "contamination"], 2.1)

    def test_validate_and_write_quality_file_accepts_names_with_extension(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            data_dir = tmp_path / "data"
            data_dir.mkdir(parents=True)
            (data_dir / "bins_to_files.json").write_text(
                json.dumps({"MAG_02": "/tmp/MAG_02.fna.gz"}),
                encoding="utf-8",
            )

            quality_file = tmp_path / "mag_qualities.tsv"
            quality_file.write_text(
                "genome\tcompleteness\tcontamination\n"
                "MAG_02.fa\t88.0\t1.5\n",
                encoding="utf-8",
            )

            self.assertTrue(cli_module.validate_and_write_quality_file(str(quality_file), str(tmp_path)))

            output = pd.read_csv(tmp_path / "cataloging" / "final" / "all_bin_metadata.csv")
            self.assertEqual(output.loc[0, "genome"], "MAG_02.fna")
            self.assertAlmostEqual(output.loc[0, "completeness"], 88.0)
            self.assertAlmostEqual(output.loc[0, "contamination"], 1.5)


if __name__ == "__main__":
    unittest.main()
