from __future__ import annotations

import builtins
import sys
from typing import Any, TextIO

try:
    from rich.console import Console
except ImportError:  # pragma: no cover - only used if runtime deps are broken.
    Console = None


console = Console(highlight=False, soft_wrap=True) if Console else None


def _target_console(file: TextIO | None):
    if Console is None:
        return None
    return Console(file=file or sys.stdout, highlight=False, soft_wrap=True)


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
    target.print(text, end=end, style=style, markup=False, highlight=False)
    if flush:
        target.file.flush()


def section(title: str) -> None:
    if console is None:
        builtins.print(f"\n{title}\n")
        return
    console.rule(f"[bold magenta]{title}[/bold magenta]", characters="-")


def prompt(message: str) -> str:
    if console is not None:
        console.print(message, end="", style="bold cyan", markup=False, highlight=False)
    else:
        builtins.print(message, end="")
    return builtins.input()
