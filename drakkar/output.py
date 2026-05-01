from __future__ import annotations

import builtins
import re
import sys
from typing import Any, TextIO

try:
    from rich.console import Console
    from rich.rule import Rule
    from rich.text import Text
    from rich.theme import Theme
except ImportError:  # pragma: no cover - only used if runtime deps are broken.
    Console = None
    Rule = None
    Text = None
    Theme = None


DRAKKAR_THEME = (
    Theme(
        {
            "drakkar.error": "bold #e85d75",
            "drakkar.info": "bold #5f9ea0",
            "drakkar.warning": "bold #d6a642",
            "drakkar.success": "bold #7fb069",
            "drakkar.heading": "bold #d6a642",
            "drakkar.rule": "#5f9ea0",
            "drakkar.prompt": "bold #5f9ea0",
            "drakkar.help": "#b7c7d3",
            "drakkar.text": "#e6edf3",
        }
    )
    if Theme
    else None
)
ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;]*m")
STATUS_STYLES = {
    "ERROR:": "drakkar.error",
    "MISSING:": "drakkar.warning",
    "INFO:": "drakkar.info",
    "PUT:": "drakkar.success",
}

console = Console(highlight=False, soft_wrap=True, theme=DRAKKAR_THEME) if Console else None


def _target_console(file: TextIO | None):
    if Console is None:
        return None
    return Console(file=file or sys.stdout, highlight=False, soft_wrap=True, theme=DRAKKAR_THEME)


def _strip_ansi(text: str) -> str:
    return ANSI_ESCAPE_RE.sub("", text)


def _status_text(text: str):
    if Text is None:
        return text

    cleaned = _strip_ansi(text)
    leading = cleaned[: len(cleaned) - len(cleaned.lstrip())]
    body = cleaned.lstrip()
    for prefix, prefix_style in STATUS_STYLES.items():
        if body.startswith(prefix):
            rendered = Text(leading)
            rendered.append(prefix, style=prefix_style)
            rendered.append(body[len(prefix) :], style="drakkar.text")
            return rendered

    stripped = cleaned.strip()
    if stripped and stripped == stripped.upper() and len(stripped.split()) <= 5:
        return Text(cleaned, style="drakkar.heading")
    if body.startswith("usage:"):
        return Text(cleaned, style="drakkar.help")
    return Text(cleaned)


def print(
    *objects: Any,
    sep: str = " ",
    end: str = "\n",
    file: TextIO | None = None,
    flush: bool = False,
    style: str | None = None,
) -> None:
    """Print through Rich while keeping the built-in print signature used here."""
    target = _target_console(file)
    if target is None:
        builtins.print(*objects, sep=sep, end=end, file=file, flush=flush)
        return

    text = sep.join(str(obj) for obj in objects)
    rendered = Text(text, style=style) if style and Text is not None else _status_text(text)
    target.print(rendered, end=end, markup=False, highlight=False)
    if flush:
        target.file.flush()


def section(title: str) -> None:
    if console is None:
        builtins.print(f"\n{title}\n")
        return
    console.print()
    console.print(Rule(Text(title, style="drakkar.heading"), characters="=", style="drakkar.rule"))


def prompt(message: str) -> str:
    if console is not None:
        console.print(message, end="", style="drakkar.prompt", markup=False, highlight=False)
    else:
        builtins.print(message, end="")
    return builtins.input()
