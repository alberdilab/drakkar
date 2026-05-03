from __future__ import annotations

import io
import re
import subprocess
import sys
import unittest
from unittest.mock import patch

from drakkar import cli as cli_module


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
        self.assertIn("Data Generation and Analysis", plain)
        self.assertIn("Operations and Management", plain)
        self.assertIn("██████████", plain)
        self.assertIn("complete", plain)
        self.assertIn("database", plain)
        self.assertIn("logging", plain)
        self.assertIn("Run the end-to-end workflow from raw reads to catalog", plain)
        self.assertIn("profiling, and annotation outputs", plain)
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

    def test_help_fallback_explains_missing_rich_dependency(self) -> None:
        parser = cli_module.RichArgumentParser(prog="drakkar-test", description="test parser")
        buffer = io.StringIO()
        with patch.object(cli_module, "rich_box", None):
            with patch.object(cli_module, "RichGroup", None):
                with patch.object(cli_module, "RichPanel", None):
                    with patch.object(cli_module, "RichTable", None):
                        with patch.object(cli_module, "RichText", None):
                            parser.print_help(buffer)

        plain = buffer.getvalue()
        self.assertIn("usage: drakkar-test", plain)
        self.assertIn("Rich-styled help is unavailable", plain)
        self.assertIn("python -m pip install rich", plain)


if __name__ == "__main__":
    unittest.main()
