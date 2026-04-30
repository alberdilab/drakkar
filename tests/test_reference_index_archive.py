from __future__ import annotations

import gzip
import importlib.util
import json
import tarfile
import tempfile
import unittest
from pathlib import Path
from subprocess import run
import sys

from drakkar.utils import file_references_to_json


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "drakkar" / "workflow" / "scripts" / "extract_reference_index.py"
BT2_SUFFIXES = [
    ".1.bt2",
    ".2.bt2",
    ".3.bt2",
    ".4.bt2",
    ".rev.1.bt2",
    ".rev.2.bt2",
]
BT2L_SUFFIXES = [
    ".1.bt2l",
    ".2.bt2l",
    ".3.bt2l",
    ".4.bt2l",
    ".rev.1.bt2l",
    ".rev.2.bt2l",
]


def load_script_module():
    spec = importlib.util.spec_from_file_location("extract_reference_index", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ReferenceIndexArchiveTests(unittest.TestCase):
    def test_extract_reference_index_renames_fasta_and_bowtie2_index(self) -> None:
        module = load_script_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            source_dir = tmp_path / "source"
            source_dir.mkdir()
            (source_dir / "host.fa").write_text(">host\nACGT\n", encoding="utf-8")
            for suffix in BT2_SUFFIXES:
                (source_dir / f"host{suffix}").write_text(f"index {suffix}\n", encoding="utf-8")

            archive_path = tmp_path / "host_index.tar.gz"
            with tarfile.open(archive_path, "w:gz") as tar:
                for path in source_dir.iterdir():
                    tar.add(path, arcname=f"nested/{path.name}")

            output_dir = tmp_path / "references"
            module.extract_reference_index(archive_path, "reference", output_dir)

            self.assertEqual((output_dir / "reference.fna").read_text(encoding="utf-8"), ">host\nACGT\n")
            for suffix in BT2_SUFFIXES:
                self.assertTrue((output_dir / f"reference{suffix}").exists())

    def test_extract_reference_index_accepts_large_bowtie2_index(self) -> None:
        module = load_script_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            source_dir = tmp_path / "source"
            source_dir.mkdir()
            (source_dir / "host.fa").write_text(">host\nACGT\n", encoding="utf-8")
            for suffix in BT2L_SUFFIXES:
                (source_dir / f"host{suffix}").write_text(f"index {suffix}\n", encoding="utf-8")

            archive_path = tmp_path / "host_index.tar.gz"
            with tarfile.open(archive_path, "w:gz") as tar:
                for path in source_dir.iterdir():
                    tar.add(path, arcname=path.name)

            output_dir = tmp_path / "references"
            module.extract_reference_index(archive_path, "reference", output_dir)

            self.assertEqual((output_dir / "reference.fna").read_text(encoding="utf-8"), ">host\nACGT\n")
            for suffix in BT2L_SUFFIXES:
                self.assertTrue((output_dir / f"reference{suffix}").exists())

    def test_extract_reference_index_decompresses_gzipped_fasta_member(self) -> None:
        module = load_script_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            source_dir = tmp_path / "source"
            source_dir.mkdir()
            with gzip.open(source_dir / "host.fna.gz", "wt", encoding="utf-8") as handle:
                handle.write(">host\nACGT\n")
            for suffix in BT2_SUFFIXES:
                (source_dir / f"host{suffix}").write_text(f"index {suffix}\n", encoding="utf-8")

            archive_path = tmp_path / "host_index.tar.gz"
            with tarfile.open(archive_path, "w:gz") as tar:
                for path in source_dir.iterdir():
                    tar.add(path, arcname=path.name)

            output_dir = tmp_path / "references"
            module.extract_reference_index(archive_path, "reference", output_dir)

            self.assertEqual((output_dir / "reference.fna").read_text(encoding="utf-8"), ">host\nACGT\n")

    def test_sample_table_reference_path_can_point_to_index_tarball(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            archive_path = tmp_path / "host_index.tar.gz"
            archive_path.write_bytes(b"placeholder")
            infofile = tmp_path / "samples.tsv"
            infofile.write_text(
                "sample\treference_name\treference_path\n"
                f"sample1\thost\t{archive_path}\n",
                encoding="utf-8",
            )

            file_references_to_json(str(infofile), tmpdir)

            reference_to_file = json.loads((tmp_path / "data" / "reference_to_file.json").read_text(encoding="utf-8"))
            self.assertEqual(reference_to_file["host"], str(archive_path.resolve()))

    def test_reference_and_reference_index_flags_are_mutually_exclusive(self) -> None:
        result = run(
            [
                sys.executable,
                "-m",
                "drakkar.cli",
                "preprocessing",
                "-r",
                "host.fna",
                "-x",
                "host_index.tar.gz",
            ],
            capture_output=True,
            text=True,
        )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not allowed with argument", result.stderr)


if __name__ == "__main__":
    unittest.main()
