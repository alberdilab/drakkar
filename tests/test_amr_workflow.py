import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar.cli_parser import build_parser
from drakkar.database_checks import amr_requirements, missing_artifacts, module_requirements
from drakkar.workflow_launcher import run_snakemake_amr


ROOT = Path(__file__).resolve().parents[1]


class AmrWorkflowTests(unittest.TestCase):
    def test_cli_exposes_amr_inputs_and_evidence_controls(self):
        args = build_parser().parse_args([
            "amr", "-f", "assemblies.tsv", "--assembly-type", "isolate",
            "--rgi-alignment-tool", "BLAST", "--rgi-include-loose",
            "--rgi-include-nudge", "--genomad-preset", "relaxed",
            "--genomad-splits", "4", "--locus-overlap", "0.75",
        ])
        self.assertEqual(args.command, "amr")
        self.assertEqual(args.assembly_type, "isolate")
        self.assertEqual(args.rgi_alignment_tool, "BLAST")
        self.assertTrue(args.rgi_include_loose)
        self.assertEqual(args.locus_overlap, 0.75)

    def test_launcher_passes_amr_configuration_to_snakemake(self):
        with patch("drakkar.workflow_launcher.run_subprocess_with_logging") as run:
            run_snakemake_amr(
                "amr", "project", "/tmp/amr-output", "/tmp/envs", "local",
                rgi_include_loose=True, genomad_splits=4, locus_overlap=0.7,
            )
        command = run.call_args.args[0][2]
        self.assertIn("workflow=amr", command)
        self.assertIn("drakkar_version=", command)
        self.assertIn("rgi_include_loose=True", command)
        self.assertIn("genomad_splits=4", command)
        self.assertIn("locus_overlap=0.7", command)

    def test_rules_call_all_three_tools_and_declare_plasmid_and_virus_outputs(self):
        rules = (ROOT / "drakkar/workflow/rules/amr.smk").read_text(encoding="utf-8")
        snakefile = (ROOT / "drakkar/workflow/Snakefile").read_text(encoding="utf-8")
        self.assertIn("amrfinder \\", rules)
        self.assertIn("rgi main \\", rules)
        self.assertIn("genomad end-to-end \\", rules)
        self.assertIn("_plasmid_summary.tsv", rules)
        self.assertIn("_virus_summary.tsv", rules)
        self.assertIn('if WORKFLOW == "amr":', snakefile)

    def test_amr_database_requirements_are_mandatory_and_card_needs_localdb(self):
        requirements = amr_requirements(config={})
        self.assertEqual(
            [requirement.config_key for requirement in requirements],
            ["AMRFINDER_DB", "CARD_DB", "GENOMAD_DB"],
        )
        self.assertIn("not configured", missing_artifacts(requirements[0])[0])

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "amrfinder").mkdir()
            (root / "amrfinder" / "version.txt").write_text("x", encoding="utf-8")
            (root / "genomad").mkdir()
            (root / "genomad" / "version.txt").write_text("x", encoding="utf-8")
            (root / "card").mkdir()
            config = {
                "AMRFINDER_DB": str(root / "amrfinder"),
                "CARD_DB": str(root / "card"),
                "GENOMAD_DB": str(root / "genomad"),
            }
            card = amr_requirements(config)[1]
            amrfinder = amr_requirements(config)[0]
            self.assertEqual(amrfinder.database, "amrfinderplus")
            self.assertEqual(card.database, "card")
            self.assertIn(str(root / "amrfinder" / "AMRProt.fa"), missing_artifacts(amrfinder))
            self.assertEqual(missing_artifacts(card), [str(root / "card" / "localDB")])
            (root / "card" / "localDB").mkdir()
            (root / "card" / "localDB" / "card.json").write_text("{}", encoding="utf-8")
            self.assertEqual(missing_artifacts(card), [])
            args = type("Args", (), {})()
            self.assertEqual(len(module_requirements("amr", args, config)), 3)

    def test_database_installer_prepares_amrfinderplus_and_card_for_runtime_use(self):
        rules = (ROOT / "drakkar/workflow/rules/databases.smk").read_text(encoding="utf-8")
        self.assertIn('if DATABASE_NAME == "amrfinderplus":', rules)
        self.assertIn("amrfinder_index", rules)
        self.assertIn('if DATABASE_NAME == "card":', rules)
        self.assertIn("rgi load --card_json", rules)
        self.assertIn("localDB/loaded_databases.json", rules)


if __name__ == "__main__":
    unittest.main()
