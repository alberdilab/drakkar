from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RULES = ROOT / "drakkar" / "workflow" / "rules" / "cataloging.smk"


def load_helpers() -> dict:
    """Exec the deep-binner resource helpers out of the Snakemake rule file."""
    text = RULES.read_text(encoding="utf-8")
    block = text.split("def _file_size_mb(path):", 1)[1].split("\nrule ", 1)[0]
    namespace: dict = {"Path": Path, "MIN_BINNING_ASSEMBLY_MB": 10}
    exec("def _file_size_mb(path):" + block, namespace)
    return namespace


def write(path: Path, size_mb: int) -> Path:
    path.write_bytes(b"0" * (size_mb * 1024 * 1024))
    return path


class DeepBinnerResourceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.helpers = load_helpers()
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def test_bam_size_drives_memory(self) -> None:
        """A small assembly with a large BAM must request far more than the floor."""
        assembly = write(self.dir / "assembly.fna", 186)
        bam = write(self.dir / "sample.bam", 1800)

        mem_mb = self.helpers["deep_binner_mem_mb"](assembly, [bam], 1)

        self.assertGreater(mem_mb, 32 * 1024)

    def test_memory_floor_and_retry_scaling(self) -> None:
        assembly = write(self.dir / "assembly.fna", 1)
        bam = write(self.dir / "sample.bam", 1)

        first = self.helpers["deep_binner_mem_mb"](assembly, [bam], 1)
        second = self.helpers["deep_binner_mem_mb"](assembly, [bam], 2)

        self.assertEqual(first, 32 * 1024)
        self.assertEqual(second, 64 * 1024)

    def test_missing_files_fall_back_to_the_floor(self) -> None:
        missing = self.dir / "absent.fna"

        self.assertEqual(self.helpers["deep_binner_mem_mb"](missing, [], 1), 32 * 1024)
        self.assertEqual(self.helpers["deep_binner_runtime"](missing, [], 1), 120)

    def test_runtime_accounts_for_every_bam(self) -> None:
        assembly = write(self.dir / "assembly.fna", 100)
        bams = [write(self.dir / f"s{i}.bam", 600) for i in range(3)]

        one = self.helpers["deep_binner_runtime"](assembly, bams[:1], 1)
        three = self.helpers["deep_binner_runtime"](assembly, bams, 1)

        self.assertGreater(three, one)


if __name__ == "__main__":
    unittest.main()
