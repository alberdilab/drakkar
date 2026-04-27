from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "drakkar" / "workflow" / "Snakefile"
CONFIG_PATH = ROOT / "drakkar" / "workflow" / "config.yaml"
RULES_DIR = ROOT / "drakkar" / "workflow" / "rules"


class WorkflowResourceTests(unittest.TestCase):
    def test_dynamic_mem_mb_resources_are_capped(self) -> None:
        uncapped_lines: list[str] = []
        pattern = re.compile(r"mem_mb\s*=\s*lambda")

        for path in sorted(RULES_DIR.glob("*.smk")):
            for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
                if pattern.search(line) and "cap_mem_mb(" not in line:
                    uncapped_lines.append(f"{path.relative_to(ROOT)}:{line_number}")

        self.assertEqual(uncapped_lines, [])

    def test_memory_cap_uses_config_with_1024_gb_default(self) -> None:
        snakefile = SNAKEFILE.read_text(encoding="utf-8")
        config = CONFIG_PATH.read_text(encoding="utf-8")
        self.assertIn('config.get("SNAKEMAKE_MAX_GB", 1024)', snakefile)
        self.assertIn("MAX_MEM_MB = SNAKEMAKE_MAX_GB * 1024", snakefile)
        self.assertIn("return min(MAX_MEM_MB, int(value))", snakefile)
        self.assertIn("SNAKEMAKE_MAX_GB: 1024", config)


if __name__ == "__main__":
    unittest.main()
