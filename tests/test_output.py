from __future__ import annotations

import contextlib
import io
import unittest

from drakkar import output


class OutputTests(unittest.TestCase):
    def test_print_strips_legacy_ansi_sequences_when_rich_is_available(self) -> None:
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            output.print("\033[1;31mERROR:\033[0m example failure")

        rendered = buffer.getvalue()
        self.assertIn("ERROR: example failure", rendered)
        self.assertNotIn("\033", rendered)

    def test_section_prints_plain_title_text(self) -> None:
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            output.section("STARTING TEST PIPELINE")

        self.assertIn("STARTING TEST PIPELINE", buffer.getvalue())


if __name__ == "__main__":
    unittest.main()
