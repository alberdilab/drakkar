from __future__ import annotations

import contextlib
import io
import sys
import unittest
from unittest.mock import patch

from drakkar import output
from drakkar import utils


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

    def test_normal_cli_streams_force_color_unless_disabled(self) -> None:
        with patch.dict(output.os.environ, {}, clear=True):
            self.assertTrue(output._should_force_terminal(sys.__stdout__))
            self.assertTrue(output._should_force_terminal(sys.__stderr__))
            self.assertFalse(output._should_force_terminal(io.StringIO()))

        with patch.dict(output.os.environ, {"NO_COLOR": "1"}, clear=True):
            self.assertTrue(output._should_force_terminal(sys.__stdout__))

        with patch.dict(output.os.environ, {"DRAKKAR_NO_COLOR": "1"}, clear=True):
            self.assertFalse(output._should_force_terminal(sys.__stdout__))

    def test_screen_session_warning_formats_command_as_code(self) -> None:
        buffer = io.StringIO()
        with (
            contextlib.redirect_stdout(buffer),
            patch.dict(utils.os.environ, {}, clear=True),
            patch("drakkar.utils.prompt", return_value="1"),
        ):
            utils.check_screen_session()

        rendered = buffer.getvalue()
        self.assertIn("To start a screen session, use:", rendered)
        self.assertIn("`screen -S mysession`", rendered)


if __name__ == "__main__":
    unittest.main()
