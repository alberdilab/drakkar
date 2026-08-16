from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from drakkar.utils import (
    check_assembly_column,
    check_reference_columns,
    detect_table_delimiter,
    file_assemblies_to_json,
    file_samples_to_json,
    read_input_table,
)


class DelimiterDetectionTests(unittest.TestCase):
    def write(self, tmpdir: str, name: str, content: str) -> str:
        path = Path(tmpdir) / name
        path.write_text(content, encoding="utf-8")
        return str(path)

    def test_detects_tabs_commas_and_semicolons(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv = self.write(tmpdir, "a.tsv", "sample\trawreads1\nsample1\tr1.fq.gz\n")
            csv = self.write(tmpdir, "a.csv", "sample,rawreads1\nsample1,r1.fq.gz\n")
            scsv = self.write(tmpdir, "a.txt", "sample;rawreads1\nsample1;r1.fq.gz\n")

            self.assertEqual(detect_table_delimiter(tsv), "\t")
            self.assertEqual(detect_table_delimiter(csv), ",")
            self.assertEqual(detect_table_delimiter(scsv), ";")

    def test_extension_does_not_decide_the_delimiter(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            mislabelled = self.write(tmpdir, "info.tsv", "sample,rawreads1\nsample1,r1.fq.gz\n")
            self.assertEqual(detect_table_delimiter(mislabelled), ",")

    def test_tabs_win_over_commas_inside_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            # The assembly column legitimately contains commas in a TSV.
            table = self.write(
                tmpdir,
                "info.tsv",
                "sample\tassembly\nsample1\tassembly1,all\n",
            )
            self.assertEqual(detect_table_delimiter(table), "\t")

    def test_falls_back_to_tab_for_single_column_and_empty_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            single = self.write(tmpdir, "single.txt", "sample\nsample1\n")
            empty = self.write(tmpdir, "empty.txt", "")

            self.assertEqual(detect_table_delimiter(single), "\t")
            self.assertEqual(detect_table_delimiter(empty), "\t")

    def test_skips_blank_leading_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            table = self.write(tmpdir, "info.csv", "\n\nsample,rawreads1\nsample1,r1.fq.gz\n")
            self.assertEqual(detect_table_delimiter(table), ",")

    def test_read_input_table_strips_column_whitespace_and_bom(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            table = self.write(tmpdir, "info.csv", "﻿sample , rawreads1\nsample1,r1.fq.gz\n")
            df = read_input_table(table, label="sample info file")
            self.assertEqual(list(df.columns), ["sample", "rawreads1"])


class CsvSampleTableTests(unittest.TestCase):
    def test_file_samples_to_json_reads_a_csv_sample_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            read1 = tmp_path / "sample1_1.fq.gz"
            read2 = tmp_path / "sample1_2.fq.gz"
            read1.write_text("r1", encoding="utf-8")
            read2.write_text("r2", encoding="utf-8")

            infofile = tmp_path / "info.csv"
            infofile.write_text(
                "sample,rawreads1,rawreads2\n"
                f"sample1,{read1},{read2}\n",
                encoding="utf-8",
            )

            file_samples_to_json(str(infofile), tmpdir)

            reads1 = json.loads((tmp_path / "data" / "sample_to_reads1.json").read_text(encoding="utf-8"))
            reads2 = json.loads((tmp_path / "data" / "sample_to_reads2.json").read_text(encoding="utf-8"))

            self.assertEqual(reads1["sample1"], [str(read1.resolve())])
            self.assertEqual(reads2["sample1"], [str(read2.resolve())])

    def test_column_probes_and_assembly_groups_work_on_a_csv_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            infofile = tmp_path / "info.csv"
            infofile.write_text(
                "sample,reference_name,reference_path,assembly\n"
                "sample1,host,host.fna,assembly1\n"
                "sample2,host,host.fna,assembly1\n",
                encoding="utf-8",
            )

            self.assertTrue(check_reference_columns(str(infofile)))
            self.assertTrue(check_assembly_column(str(infofile)))

            file_assemblies_to_json(infofile=str(infofile), output=tmpdir)
            assemblies = json.loads(
                (tmp_path / "data" / "assembly_to_samples.json").read_text(encoding="utf-8")
            )
            self.assertEqual(assemblies, {"assembly1": ["sample1", "sample2"]})

    def test_semicolon_table_is_read_as_a_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.csv"
            infofile.write_text(
                "sample;reference_name;reference_path\n"
                "sample1;host;host.fna\n",
                encoding="utf-8",
            )
            self.assertTrue(check_reference_columns(str(infofile)))


if __name__ == "__main__":
    unittest.main()
