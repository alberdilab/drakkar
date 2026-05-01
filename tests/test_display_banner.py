from __future__ import annotations

import unittest
from unittest.mock import patch

from drakkar import __version__
from drakkar import utils


class DisplayBannerTests(unittest.TestCase):
    def test_display_drakkar_embeds_version_in_ship_and_centers_intro_to_logo(self) -> None:
        with patch("drakkar.utils.print") as mocked_print:
            utils.display_drakkar()

        printed_chunks = [call.args[0] for call in mocked_print.call_args_list]
        self.assertEqual(len(printed_chunks), 3)
        ship, logo, intro = printed_chunks
        self.assertIn(f"│ v{__version__} │", ship)
        self.assertNotIn(f"│ v{__version__} │", logo)
        self.assertIn("ᚱ  Antton Alberdi  ᚱ", printed_chunks[2])
        self.assertIn("Source code: https://github.com/alberdilab/drakkar", printed_chunks[2])
        self.assertIn("╭", printed_chunks[2])
        self.assertNotIn("Version:", "\n".join(printed_chunks))
        self.assertEqual(intro, utils._intro_box(utils._ascii_block_width(logo)))
        self.assertEqual(mocked_print.call_args_list[0].kwargs["style"], utils.DRAKKAR_SHIP_STYLE)
        self.assertEqual(mocked_print.call_args_list[1].kwargs["style"], utils.DRAKKAR_LOGO_STYLE)
        self.assertEqual(mocked_print.call_args_list[2].kwargs["style"], utils.DRAKKAR_INTRO_STYLE)


if __name__ == "__main__":
    unittest.main()
