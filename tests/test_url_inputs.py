from __future__ import annotations

import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar.utils import (
    argument_references_to_json,
    file_bins_to_json,
    file_mags_to_json,
    file_references_to_json,
)


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


class UrlGenomeInputTests(unittest.TestCase):
    def test_file_references_to_json_downloads_reference_url(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            infofile = Path(tmpdir) / "info.tsv"
            infofile.write_text(
                "sample\treference_name\treference_path\n"
                "sample1\thost\thttps://example.org/reference.fna.gz\n",
                encoding="utf-8",
            )

            with patch("drakkar.utils.urlopen", return_value=FakeResponse(b">ref\nACGT\n")):
                file_references_to_json(str(infofile), tmpdir)

            output_json = Path(tmpdir) / "data" / "reference_to_file.json"
            reference_to_file = json.loads(output_json.read_text(encoding="utf-8"))
            expected_path = Path(tmpdir) / "data" / "references_cache" / "host_reference.fna.gz"

            self.assertEqual(reference_to_file["host"], str(expected_path))
            self.assertTrue(expected_path.exists())

    def test_argument_references_to_json_downloads_reference_url(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_reads = Path(tmpdir) / "sample_to_reads1.json"
            sample_reads.write_text('{"sample1": ["reads_1.fq.gz"]}', encoding="utf-8")

            with patch("drakkar.utils.urlopen", return_value=FakeResponse(b">ref\nACGT\n")):
                argument_references_to_json("https://example.org/ref.fna", str(sample_reads), tmpdir)

            output_json = Path(tmpdir) / "data" / "reference_to_file.json"
            reference_to_file = json.loads(output_json.read_text(encoding="utf-8"))
            expected_path = Path(tmpdir) / "data" / "references_cache" / "reference_ref.fna"

            self.assertEqual(reference_to_file["reference"], [str(expected_path)])
            self.assertTrue(expected_path.exists())

    def test_file_bins_to_json_downloads_genome_urls_to_cache(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            bins_file = Path(tmpdir) / "bins.txt"
            bins_file.write_text("https://example.org/bin_001.fna.gz\n", encoding="utf-8")

            with patch("drakkar.utils.urlopen", return_value=FakeResponse(b">bin\nACGT\n")):
                file_bins_to_json(str(bins_file), tmpdir)

            output_json = Path(tmpdir) / "data" / "bins_to_files.json"
            bins_to_files = json.loads(output_json.read_text(encoding="utf-8"))
            expected_path = Path(tmpdir) / "data" / "genomes_cache" / "bin_001.fna.gz"

            self.assertEqual(bins_to_files["bin_001"], str(expected_path))
            self.assertTrue(expected_path.exists())

    def test_file_mags_to_json_downloads_genome_urls_and_normalizes_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            mags_file = Path(tmpdir) / "mags.txt"
            mags_file.write_text("https://example.org/MAG_01.fna.gz\n", encoding="utf-8")

            with patch("drakkar.utils.urlopen", return_value=FakeResponse(b">mag\nACGT\n")):
                file_mags_to_json(str(mags_file), tmpdir)

            output_json = Path(tmpdir) / "data" / "mags_to_files.json"
            mags_to_files = json.loads(output_json.read_text(encoding="utf-8"))
            expected_path = Path(tmpdir) / "data" / "genomes_cache" / "MAG_01.fna.gz"

            self.assertEqual(mags_to_files["MAG_01"], str(expected_path))
            self.assertTrue(expected_path.exists())


if __name__ == "__main__":
    unittest.main()
