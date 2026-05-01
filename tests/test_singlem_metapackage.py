from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RULES_DIR = ROOT / "drakkar" / "workflow" / "rules"
SINGLEM_COMMANDS = (
    "singlem pipe",
    "singlem microbial_fraction",
    "singlem prokaryotic_fraction",
)


class SinglemMetapackageTests(unittest.TestCase):
    def test_all_singlem_workflow_rules_use_configured_metapackage(self) -> None:
        missing_metapackage: list[str] = []
        singlem_rule_count = 0

        for path in sorted(RULES_DIR.glob("*.smk")):
            blocks = re.split(r"(?m)^rule\s+", path.read_text(encoding="utf-8"))
            for block in blocks:
                if any(command in block for command in SINGLEM_COMMANDS):
                    singlem_rule_count += 1
                    if "--metapackage {params.singlem_db}" not in block:
                        rule_name = block.split(":", maxsplit=1)[0].strip()
                        missing_metapackage.append(f"{path.relative_to(ROOT)}:{rule_name}")

        self.assertEqual(missing_metapackage, [])
        self.assertGreaterEqual(singlem_rule_count, 4)


if __name__ == "__main__":
    unittest.main()
