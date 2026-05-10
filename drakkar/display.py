import os
import sys
import time

from drakkar import __version__
from drakkar.ascii import (
    DRAKKAR_LOGO_ART,
    DRAKKAR_SHIP_ART,
    END_ART,
    INTRO_ROWS,
    UNLOCK_ART,
    UPDATE_SUCCESS_ASCII_TEMPLATE,
)
from drakkar.output import Text as RichText, print

DRAKKAR_SHIP_STYLE = "bold #5f9ea0"

DRAKKAR_LOGO_STYLE = "bold #d6a642"

DRAKKAR_INTRO_STYLE = "bold #b7c7d3"

DRAKKAR_VERSION_BADGE_STYLE = DRAKKAR_INTRO_STYLE

BANNER_DELAY_SECONDS = 0.3

def _version_badge_lines(version):
    label = f"v{version}"
    width = len(label) + 2
    return [
        "╭" + "─" * width + "╮",
        f"│ {label} │",
        "╰" + "─" * width + "╯",
    ]

def _ascii_block_with_bottom_right_badge(ascii_text, version):
    lines = ascii_text.splitlines()
    if not lines:
        return ascii_text

    non_empty_line_indexes = [index for index, line in enumerate(lines) if line.strip()]
    if not non_empty_line_indexes:
        return ascii_text

    badge_lines = _version_badge_lines(version)
    badge_width = max(len(line) for line in badge_lines)
    block_width = max(len(line) for line in lines)
    badge_start = max(0, block_width - badge_width)
    first_badge_line = max(
        non_empty_line_indexes[0],
        non_empty_line_indexes[-1] - len(badge_lines) + 1,
    )
    target_indexes = [first_badge_line + offset for offset in range(len(badge_lines))]
    content_width = max(
        len(lines[index].rstrip()) for index in target_indexes if index < len(lines)
    )
    if content_width > badge_start:
        badge_start = content_width + 2

    for offset, badge_line in enumerate(badge_lines):
        line_index = first_badge_line + offset
        if line_index >= len(lines):
            lines.append("")
        lines[line_index] = lines[line_index].ljust(badge_start) + badge_line

    return "\n".join(lines)

def _ascii_block_width(ascii_text):
    return max((len(line) for line in ascii_text.splitlines()), default=0)

def _center_block(block, target_width):
    block_width = _ascii_block_width(block)
    padding = max(0, (target_width - block_width) // 2)
    return "\n".join(
        (" " * padding + line) if line.strip() else line
        for line in block.splitlines()
    )

def _intro_box(target_width):
    content_width = max(len(row) for row in INTRO_ROWS)
    inner_width = content_width + 4
    lines = [
        "╭" + "─" * inner_width + "╮",
        *[f"│  {row.center(content_width)}  │" for row in INTRO_ROWS],
        "╰" + "─" * inner_width + "╯",
    ]
    return _center_block("\n".join(lines), target_width)

def _styled_ascii_art(ascii_text, base_style, accent_chunks=None, accent_style=None):
    if RichText is None:
        return ascii_text

    rendered = RichText(ascii_text, style=base_style, no_wrap=True, overflow="ignore")
    if accent_chunks and accent_style:
        search_start = 0
        for chunk in accent_chunks:
            start = ascii_text.find(chunk, search_start)
            if start == -1:
                start = ascii_text.find(chunk)
            if start == -1:
                continue
            rendered.stylize(accent_style, start, start + len(chunk))
            search_start = start + len(chunk)
    return rendered

def get_drakkar_banner_blocks(include_intro=True):
    ascii_ship = _ascii_block_with_bottom_right_badge(DRAKKAR_SHIP_ART, __version__)
    blocks = [
        (ascii_ship, DRAKKAR_SHIP_STYLE),
        (DRAKKAR_LOGO_ART, DRAKKAR_LOGO_STYLE),
    ]
    if include_intro:
        ascii_intro = _intro_box(_ascii_block_width(DRAKKAR_LOGO_ART))
        blocks.append((ascii_intro, DRAKKAR_INTRO_STYLE))
    return blocks

def get_drakkar_banner_renderables(include_intro=True):
    ascii_ship = _ascii_block_with_bottom_right_badge(DRAKKAR_SHIP_ART, __version__)
    renderables = [
        _styled_ascii_art(
            ascii_ship,
            DRAKKAR_SHIP_STYLE,
            accent_chunks=_version_badge_lines(__version__),
            accent_style=DRAKKAR_VERSION_BADGE_STYLE,
        ),
        _styled_ascii_art(DRAKKAR_LOGO_ART, DRAKKAR_LOGO_STYLE),
    ]
    if include_intro:
        ascii_intro = _intro_box(_ascii_block_width(DRAKKAR_LOGO_ART))
        renderables.append(_styled_ascii_art(ascii_intro, DRAKKAR_INTRO_STYLE))
    return renderables

def _banner_animation_enabled(stream=None):
    if os.environ.get("DRAKKAR_NO_ANIMATION"):
        return False
    stream = stream or sys.stdout
    return bool(getattr(stream, "isatty", lambda: False)())

def display_banner_sequence(renderables, file=None, delay_after=False):
    delay_enabled = _banner_animation_enabled(file)
    rendered = list(renderables)
    for index, renderable in enumerate(rendered):
        is_last = index == len(rendered) - 1
        print(renderable, file=file)
        if delay_enabled and (not is_last or delay_after):
            time.sleep(BANNER_DELAY_SECONDS)

def display_drakkar(file=None):
    display_banner_sequence(get_drakkar_banner_renderables(), file=file)

def display_unlock():
    print(UNLOCK_ART)

def display_end():
    print(END_ART)

def _format_update_success_ascii(version):
    return UPDATE_SUCCESS_ASCII_TEMPLATE.replace("X.X.XX", str(version).center(len("X.X.XX")), 1)

def display_update_success(version):
    print(_styled_ascii_art(_format_update_success_ascii(version), DRAKKAR_SHIP_STYLE))
