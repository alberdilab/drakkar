from __future__ import annotations

import contextlib
import io
import re
import unittest
from pathlib import Path

from drakkar.cli import normalize_annotation_type
from drakkar.cli_parser import build_parser


ROOT = Path(__file__).resolve().parents[1]
ANNOTATION_RULES = ROOT / "drakkar" / "workflow" / "rules" / "annotating_function.smk"
DBCAN_ENV = ROOT / "drakkar" / "workflow" / "envs" / "annotating_function_dbcan.yaml"


class AnnotationTableWorkflowTests(unittest.TestCase):
    def test_merge_gene_annotations_passes_filter_params(self) -> None:
        rules = ANNOTATION_RULES.read_text(encoding="utf-8")

        self.assertIn('ANNOTATION_EVALUE = float(config.get("annotation_evalue"', rules)
        self.assertIn('ANNOTATION_IDENTITY = float(config.get("annotation_identity"', rules)
        self.assertIn('config.get("ANNOTATION_IDENTITY", 50.0)', rules)
        self.assertIn('"annotation_query_coverage", config.get("ANNOTATION_QUERY_COVERAGE", 0.5)', rules)
        self.assertIn('"annotation_target_coverage", config.get("ANNOTATION_TARGET_COVERAGE", 0.5)', rules)
        self.assertIn("--evalue {params.evalue}", rules)
        self.assertIn("--identity {params.identity}", rules)
        self.assertIn("--query-coverage {params.query_coverage}", rules)
        self.assertIn("--target-coverage {params.target_coverage}", rules)
        self.assertIn("--sources {params.sources}", rules)
        self.assertIn("--qc-output {output.qc}", rules)
        self.assertIn("--mag {params.mag}", rules)

    def test_vfdb_rule_emits_query_and_target_coverage(self) -> None:
        rules = ANNOTATION_RULES.read_text(encoding="utf-8")
        match = re.search(r"rule vfdb:.*?(?=\nrule )", rules, re.DOTALL)

        self.assertIsNotNone(match)
        vfdb_rule = match.group(0)
        self.assertIn("qlen,tlen,qcov,tcov", vfdb_rule)

    def test_annotation_report_is_a_final_functional_output(self) -> None:
        rules = ANNOTATION_RULES.read_text(encoding="utf-8")
        snakefile = (ROOT / "drakkar" / "workflow" / "Snakefile").read_text(encoding="utf-8")

        self.assertIn("rule annotation_report:", rules)
        self.assertIn("annotation_manifest.yaml", snakefile)
        self.assertIn("annotation_qc.tsv", snakefile)

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

    def test_cazy_rule_uses_dbcan_coverage_filtered_output(self) -> None:
        rules = ANNOTATION_RULES.read_text(encoding="utf-8")
        match = re.search(r"rule cazy:.*?(?=\nrule )", rules, re.DOTALL)

        self.assertIsNotNone(match)
        cazy_rule = match.group(0)
        self.assertIn("dbCAN_hmm_results.tsv", cazy_rule)
        self.assertIn("run_dbcan CAZyme_annotation", cazy_rule)
        self.assertIn("--mode protein", cazy_rule)
        self.assertIn("--methods hmm", cazy_rule)
        self.assertIn("evalue=1e-15", cazy_rule)
        self.assertIn("coverage=0.35", cazy_rule)
        self.assertIn("--e_value_threshold_dbcan {params.evalue}", cazy_rule)
        self.assertIn("--coverage_threshold_dbcan {params.coverage}", cazy_rule)
        self.assertIn("if [ ! -s \"{output}\" ]", cazy_rule)
        self.assertNotIn("hmmscan", cazy_rule)

        self.assertIn(
            'selected.append(f"{OUTPUT_DIR}/annotating/cazy/{wildcards.mag}/dbCAN_hmm_results.tsv")',
            rules,
        )

    def test_kegg_rule_leaves_acceptance_to_native_kofam_cutoffs(self) -> None:
        rules = ANNOTATION_RULES.read_text(encoding="utf-8")
        match = re.search(r"rule kegg:.*?(?=\nrule )", rules, re.DOTALL)

        self.assertIsNotNone(match)
        kegg_rule = match.group(0)
        self.assertIn("-E 10 --domE 10", kegg_rule)
        self.assertNotIn("-E 1e-10", kegg_rule)

    def test_structure_annotation_is_work_in_progress_and_unavailable(self) -> None:
        for requested in ("structure", "foldseek", "kegg,structure"):
            with self.subTest(requested=requested):
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    normalized = normalize_annotation_type(requested)
                self.assertIsNone(normalized)
                self.assertIn("work in progress", output.getvalue())
                self.assertIn("not available in Drakkar 2.0.0", output.getvalue())

        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                build_parser().parse_args(["database", "foldseek"])

    def test_dbcan_environment_pins_compatible_pyhmmer_api(self) -> None:
        environment = DBCAN_ENV.read_text(encoding="utf-8")

        self.assertIn("- python=3.12", environment)
        self.assertIn("- dbcan=5.2.5", environment)
        self.assertIn("- pyhmmer=0.11.4", environment)


if __name__ == "__main__":
    unittest.main()
