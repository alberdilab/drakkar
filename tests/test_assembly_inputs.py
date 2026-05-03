from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from drakkar.utils import check_assembly_column, file_assemblies_to_json


class AssemblyInputTests(unittest.TestCase):
    def test_check_assembly_column_accepts_preferred_assembly_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            table = Path(tmpdir) / "samples.tsv"
            table.write_text(
                "sample\tassembly\n"
                "sample1\tassembly_alpha\n"
                "sample2\tassembly_beta\n",
                encoding="utf-8",
            )

            self.assertTrue(check_assembly_column(str(table)))

    def test_check_assembly_column_accepts_legacy_coassembly_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            table = Path(tmpdir) / "samples.tsv"
            table.write_text(
                "sample\tcoassembly\n"
                "sample1\tassembly_alpha\n",
                encoding="utf-8",
            )

            self.assertTrue(check_assembly_column(str(table)))

    def test_file_assemblies_to_json_prefers_assembly_column_over_legacy_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            table = Path(tmpdir) / "samples.tsv"
            table.write_text(
                "sample\tassembly\tcoassembly\n"
                "sample1\tassembly_alpha\tlegacy_one\n"
                "sample2\tassembly_beta\tlegacy_two\n",
                encoding="utf-8",
            )

            file_assemblies_to_json(str(table), output=tmpdir)

            with open(Path(tmpdir) / "data" / "assembly_to_samples.json", "r", encoding="utf-8") as handle:
                mapping = json.load(handle)

            self.assertEqual(
                mapping,
                {
                    "assembly_alpha": ["sample1"],
                    "assembly_beta": ["sample2"],
                },
            )

    def test_file_assemblies_to_json_uses_distinct_assembly_values_as_explicit_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            table = Path(tmpdir) / "samples.tsv"
            table.write_text(
                "sample\tassembly\n"
                "sampleA\tassembly_alpha\n"
                "sampleB\tassembly_beta\n",
                encoding="utf-8",
            )

            file_assemblies_to_json(str(table), output=tmpdir)

            with open(Path(tmpdir) / "data" / "assembly_to_samples.json", "r", encoding="utf-8") as handle:
                mapping = json.load(handle)

            self.assertEqual(
                mapping,
                {
                    "assembly_alpha": ["sampleA"],
                    "assembly_beta": ["sampleB"],
                },
            )

    def test_file_assemblies_to_json_adds_individual_mode_alongside_named_assemblies(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            table = Path(tmpdir) / "samples.tsv"
            table.write_text(
                "sample\tassembly\n"
                "sampleA\tassembly_alpha\n"
                "sampleB\tassembly_beta\n",
                encoding="utf-8",
            )

            file_assemblies_to_json(
                str(table),
                samples=["sampleA", "sampleB"],
                individual=True,
                output=tmpdir,
            )

            with open(Path(tmpdir) / "data" / "assembly_to_samples.json", "r", encoding="utf-8") as handle:
                mapping = json.load(handle)

            self.assertEqual(
                mapping,
                {
                    "assembly_alpha": ["sampleA"],
                    "assembly_beta": ["sampleB"],
                    "sampleA": ["sampleA"],
                    "sampleB": ["sampleB"],
                },
            )


if __name__ == "__main__":
    unittest.main()
