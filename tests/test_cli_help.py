from __future__ import annotations

import re
import subprocess
import sys
import unittest


ANSI_RE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


class CliHelpTests(unittest.TestCase):
    def run_help(self, *args: str) -> str:
        result = subprocess.run(
            [sys.executable, "-m", "drakkar.cli", *args, "--help"],
            capture_output=True,
            text=True,
            check=True,
        )
        return ANSI_RE.sub("", result.stdout)

    def test_top_level_help_groups_commands(self) -> None:
        plain = self.run_help()
        self.assertIn("Workflow launcher and operations hub", plain)
        self.assertIn("Start Here", plain)
        self.assertIn("Data Generation Workflows", plain)
        self.assertIn("Analysis Workflows", plain)
        self.assertIn("Operations and Management", plain)
        self.assertIn("complete", plain)
        self.assertIn("database", plain)
        self.assertIn("logging", plain)
        self.assertIn("Examples", plain)
        self.assertNotIn("\noptions:\n", plain.lower())

    def test_preprocessing_help_uses_grouped_sections(self) -> None:
        plain = self.run_help("preprocessing")
        self.assertIn("drakkar preprocessing", plain)
        self.assertIn("Category", plain)
        self.assertIn("Examples", plain)
        self.assertIn("Input Sources", plain)
        self.assertIn("Optional Analyses", plain)
        self.assertIn("Run Configuration", plain)
        self.assertIn("Resource Scaling", plain)
        self.assertIn("Option", plain)
        self.assertIn("Value", plain)
        self.assertIn("Description", plain)
        self.assertIn("-i, --input", plain)
        self.assertIn("--fraction", plain)
        self.assertNotIn("\noptions:\n", plain.lower())

    def test_profiling_help_uses_grouped_sections(self) -> None:
        plain = self.run_help("profiling")
        self.assertIn("Analysis workflow", plain)
        self.assertIn("Input Sources", plain)
        self.assertIn("Analysis Settings", plain)
        self.assertIn("Run Configuration", plain)
        self.assertIn("Resource Scaling", plain)
        self.assertIn("-b, --bins_dir", plain)
        self.assertIn("-n, --ignore_quality", plain)
        self.assertIn("--memory-multiplier", plain)
        self.assertIn("Examples", plain)


if __name__ == "__main__":
    unittest.main()
