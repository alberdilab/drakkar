from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "drakkar" / "workflow" / "Snakefile"
CONFIG_PATH = ROOT / "drakkar" / "workflow" / "config.yaml"
RULES_DIR = ROOT / "drakkar" / "workflow" / "rules"


class WorkflowResourceTests(unittest.TestCase):
    def test_mem_mb_resources_are_capped(self) -> None:
        uncapped_lines: list[str] = []
        pattern = re.compile(r"mem_mb\s*=")

        for path in sorted(RULES_DIR.glob("*.smk")):
            for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
                if pattern.search(line) and "cap_mem_mb(" not in line:
                    uncapped_lines.append(f"{path.relative_to(ROOT)}:{line_number}")

        self.assertEqual(uncapped_lines, [])

    def test_runtime_resources_are_capped(self) -> None:
        uncapped_lines: list[str] = []
        pattern = re.compile(r"runtime\s*=")

        for path in sorted(RULES_DIR.glob("*.smk")):
            for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
                if pattern.search(line) and "cap_runtime(" not in line:
                    uncapped_lines.append(f"{path.relative_to(ROOT)}:{line_number}")

        self.assertEqual(uncapped_lines, [])

    def test_resource_caps_use_config_defaults_and_multipliers(self) -> None:
        snakefile = SNAKEFILE.read_text(encoding="utf-8")
        config = CONFIG_PATH.read_text(encoding="utf-8")
        self.assertIn('positive_config_int("SNAKEMAKE_MAX_GB", 1024)', snakefile)
        self.assertIn('positive_config_int("SNAKEMAKE_MAX_TIME", 14 * 24 * 60)', snakefile)
        self.assertIn('config.get("MEMORY_MULTIPLIER", 1)', snakefile)
        self.assertIn('config.get("TIME_MULTIPLIER", 1)', snakefile)
        self.assertIn("MAX_MEM_MB = SNAKEMAKE_MAX_GB * 1024", snakefile)
        self.assertIn("MAX_RUNTIME_MIN = SNAKEMAKE_MAX_TIME", snakefile)
        self.assertIn("return min(MAX_MEM_MB, int(math.ceil(float(value))) * MEMORY_MULTIPLIER)", snakefile)
        self.assertIn("return min(MAX_RUNTIME_MIN, int(math.ceil(float(value))) * TIME_MULTIPLIER)", snakefile)
        self.assertIn("SNAKEMAKE_MAX_GB: 1024", config)
        self.assertIn("SNAKEMAKE_MAX_TIME: 20160", config)
        self.assertIn("MEMORY_MULTIPLIER: 1", config)
        self.assertIn("TIME_MULTIPLIER: 1", config)


if __name__ == "__main__":
    unittest.main()
