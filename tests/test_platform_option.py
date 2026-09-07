from __future__ import annotations

import re
import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar import cli as cli_module
from drakkar.cli_parser import build_parser

RULES_DIR = Path(__file__).resolve().parents[1] / "drakkar" / "workflow" / "rules"

ILLUMINA_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ILLUMINA_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
BGI_R1 = "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"
BGI_R2 = "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"


def platform_header(rule_file: str, marker: str, platform: str) -> dict:
    """Execute the platform block at the top of a rules file with a fake config."""
    source = (RULES_DIR / rule_file).read_text()
    block = source.split(marker)[1].split("####\n# Workflow rules\n####")[0]
    namespace = {"config": {"platform": platform}}
    exec(block, namespace)
    return namespace


class PlatformArgumentTests(unittest.TestCase):
    def test_platform_defaults_to_illumina(self) -> None:
        for command in ("preprocessing", "complete"):
            with self.subTest(command=command):
                args = build_parser().parse_args([command, "-f", "samples.tsv"])
                self.assertEqual(args.platform, "illumina")

    def test_platform_accepts_bgi(self) -> None:
        for command in ("preprocessing", "complete"):
            with self.subTest(command=command):
                args = build_parser().parse_args([command, "-f", "samples.tsv", "--platform", "bgi"])
                self.assertEqual(args.platform, "bgi")

    def test_unknown_platform_is_rejected(self) -> None:
        with self.assertRaises(SystemExit):
            build_parser().parse_args(["preprocessing", "-f", "samples.tsv", "--platform", "nanopore"])


class PlatformPropagationTests(unittest.TestCase):
    def command_for_platform(self, platform: str) -> str:
        with patch.object(cli_module, "config_vars", {"SNAKEMAKE_MODULE": "snakemake"}):
            with patch.object(cli_module, "run_subprocess_with_logging") as run_mock:
                cli_module.run_snakemake_preprocessing(
                    "preprocessing", "project", "/tmp/output", False, "/tmp/envs", "local",
                    False, False, False, platform, 1, 1,
                )
        return run_mock.call_args.args[0][2]

    def test_platform_reaches_snakemake_config(self) -> None:
        for platform in ("illumina", "bgi"):
            with self.subTest(platform=platform):
                self.assertIn(f"platform={platform}", self.command_for_platform(platform))


class PlatformAdapterTests(unittest.TestCase):
    RULE_FILES = ("preprocessing.smk", "preprocessing_ref.smk")
    MARKER = "####\n# Platform-dependent read handling\n####\n"

    def test_adapters_match_platform(self) -> None:
        expected = {
            "illumina": (ILLUMINA_R1, ILLUMINA_R2),
            "bgi": (BGI_R1, BGI_R2),
        }
        for rule_file in self.RULE_FILES:
            for platform, (r1, r2) in expected.items():
                with self.subTest(rule_file=rule_file, platform=platform):
                    ns = platform_header(rule_file, self.MARKER, platform)
                    self.assertEqual(ns["ADAPTER_R1"], r1)
                    self.assertEqual(ns["ADAPTER_R2"], r2)

    def test_polyg_trimming_only_for_illumina(self) -> None:
        for rule_file in self.RULE_FILES:
            with self.subTest(rule_file=rule_file):
                self.assertEqual(
                    platform_header(rule_file, self.MARKER, "illumina")["POLYG_FLAG"], "--trim_poly_g"
                )
                self.assertEqual(platform_header(rule_file, self.MARKER, "bgi")["POLYG_FLAG"], "")

    def test_unknown_platform_raises(self) -> None:
        for rule_file in self.RULE_FILES:
            with self.subTest(rule_file=rule_file):
                with self.assertRaises(ValueError):
                    platform_header(rule_file, self.MARKER, "nanopore")

    def test_no_hardcoded_adapters_remain_in_shell(self) -> None:
        for rule_file in self.RULE_FILES:
            with self.subTest(rule_file=rule_file):
                shell = re.search(
                    r'rule fastp:.*?shell:\n        """\n(.*?)\n        """',
                    (RULES_DIR / rule_file).read_text(),
                    re.S,
                ).group(1)
                self.assertNotIn(ILLUMINA_R1, shell)
                self.assertNotIn(ILLUMINA_R2, shell)
                self.assertIn("{params.adapter_r1}", shell)
                self.assertIn("{params.adapter_r2}", shell)
                self.assertIn("--detect_adapter_for_pe", shell)


class SanitizeReadNameTests(unittest.TestCase):
    MARKER = 'SEQKIT_MODULE = config["SEQKIT_MODULE"]\n'

    def test_id_regexp_only_for_bgi(self) -> None:
        self.assertEqual(
            platform_header("preprocessing_sanitize.smk", self.MARKER, "illumina")["PAIR_ID_REGEXP"], ""
        )
        self.assertIn(
            "--id-regexp",
            platform_header("preprocessing_sanitize.smk", self.MARKER, "bgi")["PAIR_ID_REGEXP"],
        )

    def test_id_regexp_pairs_bgi_mates(self) -> None:
        """The regexp must strip the /1 and /2 suffix so mates share an ID."""
        flag = platform_header("preprocessing_sanitize.smk", self.MARKER, "bgi")["PAIR_ID_REGEXP"]
        pattern = re.compile(flag.split("'")[1])
        for stem in ("V300026712L2C001R0010000372", "CL100050407L1C001R001_1"):
            with self.subTest(stem=stem):
                self.assertEqual(pattern.match(f"{stem}/1").group(1), stem)
                self.assertEqual(pattern.match(f"{stem}/2").group(1), stem)

    def test_zero_pair_guard_present(self) -> None:
        source = (RULES_DIR / "preprocessing_sanitize.smk").read_text()
        self.assertIn('if [ "$count" -eq 0 ]; then', source)
        self.assertIn("exit 1", source)


if __name__ == "__main__":
    unittest.main()
