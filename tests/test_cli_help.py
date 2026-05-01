from __future__ import annotations

import re
import subprocess
import sys
import unittest


ANSI_RE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


class CliHelpTests(unittest.TestCase):
    def test_preprocessing_help_uses_rich_tables(self) -> None:
        result = subprocess.run(
            [sys.executable, "-m", "drakkar.cli", "preprocessing", "--help"],
            capture_output=True,
            text=True,
            check=True,
        )

        plain = ANSI_RE.sub("", result.stdout)
        self.assertIn("drakkar preprocessing", plain)
        self.assertIn("Options", plain)
        self.assertIn("Option", plain)
        self.assertIn("Value", plain)
        self.assertIn("Description", plain)
        self.assertIn("-i, --input", plain)
        self.assertIn("--fraction", plain)
        self.assertNotIn("\noptions:\n", plain.lower())


if __name__ == "__main__":
    unittest.main()
