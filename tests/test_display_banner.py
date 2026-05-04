from __future__ import annotations

import unittest
from unittest.mock import patch

from drakkar import __version__
from drakkar import utils


class DisplayBannerTests(unittest.TestCase):
    def test_display_drakkar_embeds_version_in_ship_and_centers_intro_to_logo(self) -> None:
        with patch("drakkar.utils.print") as mocked_print:
            utils.display_drakkar()

        printed_chunks = [getattr(call.args[0], "plain", call.args[0]) for call in mocked_print.call_args_list]
        self.assertEqual(len(printed_chunks), 3)
        ship, logo, intro = printed_chunks
        self.assertIn(f"│ v{__version__} │", ship)
        self.assertNotIn(f"│ v{__version__} │", logo)
        self.assertIn("ᚱ  Antton Alberdi  ᚱ", printed_chunks[2])
        self.assertIn("Source code: https://github.com/alberdilab/drakkar", printed_chunks[2])
        self.assertIn("╭", printed_chunks[2])
        self.assertNotIn("Version:", "\n".join(printed_chunks))
        self.assertEqual(intro, utils._intro_box(utils._ascii_block_width(logo)))
        ship_renderable = mocked_print.call_args_list[0].args[0]
        self.assertTrue(any(span.style == utils.DRAKKAR_VERSION_BADGE_STYLE for span in ship_renderable.spans))

    def test_display_update_success_replaces_version_placeholder(self) -> None:
        with patch("drakkar.utils.print") as mocked_print:
            utils.display_update_success("1.4.0")

        banner = getattr(mocked_print.call_args.args[0], "plain", mocked_print.call_args.args[0])
        self.assertIn("You have succesfully installed", banner)
        self.assertIn("DRAKKAR version", banner)
        self.assertIn("1.4.0", banner)
        self.assertNotIn("X.X.XX", banner)

    def test_display_drakkar_pauses_between_banner_blocks_when_animation_enabled(self) -> None:
        with (
            patch("drakkar.utils.print"),
            patch("drakkar.utils._banner_animation_enabled", return_value=True),
            patch("drakkar.utils.time.sleep") as mocked_sleep,
        ):
            utils.display_drakkar()

        self.assertEqual(mocked_sleep.call_count, 2)
        mocked_sleep.assert_called_with(utils.BANNER_DELAY_SECONDS)


if __name__ == "__main__":
    unittest.main()
