from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANNOTATION_RULES = ROOT / "drakkar" / "workflow" / "rules" / "annotating_function.smk"


class AnnotationTableWorkflowTests(unittest.TestCase):
    def test_final_annotation_tables_strip_duplicate_per_genome_headers(self) -> None:
        rules = ANNOTATION_RULES.read_text(encoding="utf-8")

        for rule_name in ("final_gene_annotation_table", "final_cluster_annotation_table"):
            with self.subTest(rule=rule_name):
                match = re.search(
                    rf"rule {rule_name}:.*?shell:\s*\"\"\"(?P<shell>.*?)\"\"\"",
                    rules,
                    re.DOTALL,
                )
                self.assertIsNotNone(match)
                shell = match.group("shell")
                self.assertIn("awk 'FNR==1 && NR!=1 {{ next }} {{ print }}'", shell)
                self.assertIn("| xz -c > {output}", shell)


if __name__ == "__main__":
    unittest.main()
