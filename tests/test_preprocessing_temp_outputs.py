from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PREPROCESSING_REF_RULES = ROOT / "drakkar" / "workflow" / "rules" / "preprocessing_ref.smk"


class PreprocessingTempOutputsTests(unittest.TestCase):
    def test_host_and_metagenomic_count_sidecars_are_temporary(self) -> None:
        rules = PREPROCESSING_REF_RULES.read_text(encoding="utf-8")
        expected_outputs = [
            "metareads",
            "metabases",
            "hostreads",
            "hostbases",
        ]

        for output_name in expected_outputs:
            self.assertIn(
                f'{output_name}=temp(f"{{OUTPUT_DIR}}/preprocessing/final/{{{{sample}}}}.{output_name}")',
                rules,
            )

    def test_preprocessing_stats_still_consumes_count_sidecars(self) -> None:
        rules = PREPROCESSING_REF_RULES.read_text(encoding="utf-8")
        stats_rule = re.search(r"rule preprocessing_stats:.*?(?=\nrule |\Z)", rules, re.S)

        self.assertIsNotNone(stats_rule)
        stats_block = stats_rule.group(0)
        for suffix in [".metareads", ".metabases", ".hostreads", ".hostbases"]:
            self.assertIn(suffix, stats_block)


if __name__ == "__main__":
    unittest.main()
