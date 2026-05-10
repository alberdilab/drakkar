import argparse
import sys

try:
    from rich import box as rich_box
    from rich.console import Group as RichGroup
    from rich.panel import Panel as RichPanel
    from rich.table import Table as RichTable
    from rich.text import Text as RichText
except ImportError:  # pragma: no cover - only used if runtime deps are broken.
    rich_box = None
    RichGroup = None
    RichPanel = None
    RichTable = None
    RichText = None

from drakkar.display import display_banner_sequence, get_drakkar_banner_renderables
from drakkar.output import get_console, print

def _rich_available():
    return all((rich_box, RichGroup, RichPanel, RichTable, RichText))

def _clean_usage(parser):
    try:
        formatter = parser._get_formatter()
        formatter._width = 68
        formatter.add_usage(
            parser.usage,
            parser._actions,
            parser._mutually_exclusive_groups,
        )
        usage = formatter.format_help().strip()
    except Exception:
        usage = parser.format_usage().strip()
    if usage.startswith("usage: "):
        usage = usage[len("usage: "):]
    return usage

def _action_value(action):
    if getattr(action, "nargs", None) == 0:
        return ""
    if action.metavar is not None:
        if isinstance(action.metavar, tuple):
            value = " ".join(str(item) for item in action.metavar)
        else:
            value = str(action.metavar)
    elif action.choices is not None:
        value = "{" + ",".join(str(choice) for choice in action.choices) + "}"
    else:
        value = str(action.dest).upper()

    if action.nargs == "+":
        return f"{value} ..."
    if action.nargs == "*":
        return f"[{value} ...]"
    if action.nargs == "?":
        return f"[{value}]"
    return value

def _expanded_action_help(parser, action):
    if action.help in (None, argparse.SUPPRESS):
        return ""
    try:
        help_text = parser._get_formatter()._expand_help(action)
    except Exception:
        help_text = str(action.help)
    if getattr(action, "required", False):
        help_text = f"{help_text} Required." if help_text else "Required."
    return help_text

def _display_group_title(title):
    if title in {"options", "optional arguments"}:
        return "Options"
    if title == "positional arguments":
        return "Arguments"
    return str(title).replace("_", " ").title()

def _set_help_metadata(parser, *, category=None, examples=None, sections=None, command_groups=None):
    parser._help_category = category
    parser._help_examples = list(examples or [])
    parser._help_sections = list(sections or [])
    parser._help_command_groups = list(command_groups or [])
    return parser

def _rich_table(title):
    table = RichTable(
        title=RichText(title, style="drakkar.heading"),
        box=rich_box.SIMPLE_HEAVY,
        border_style="drakkar.rule",
        header_style="drakkar.heading",
        expand=True,
        padding=(0, 1),
    )
    table.add_column("Option", style="bold #d6a642", no_wrap=True)
    table.add_column("Value", style="#5f9ea0", no_wrap=True)
    table.add_column("Description", style="drakkar.help")
    return table

def _find_action_by_key(parser, key):
    for action in parser._actions:
        if action.dest == key:
            return action
        if key in getattr(action, "option_strings", ()):
            return action
    return None

def _rich_action_table(parser, action_group):
    actions = [
        action
        for action in action_group._group_actions
        if action.help is not argparse.SUPPRESS
        and not isinstance(action, argparse._SubParsersAction)
    ]
    if not actions:
        return None

    table = _rich_table(_display_group_title(action_group.title))
    for action in actions:
        option = ", ".join(action.option_strings) if action.option_strings else action.dest
        table.add_row(option, _action_value(action), _expanded_action_help(parser, action))
    return table

def _subparser_action(parser):
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action
    return None

def _rich_command_table(title, rows):
    if not rows:
        return None
    table = RichTable(
        title=RichText(title, style="drakkar.heading"),
        box=rich_box.SIMPLE_HEAVY,
        border_style="drakkar.rule",
        header_style="drakkar.heading",
        expand=True,
        padding=(0, 1),
    )
    table.add_column("Command", style="bold #d6a642", no_wrap=True)
    table.add_column("Purpose", style="drakkar.help")
    for command, purpose in rows:
        table.add_row(command, purpose or "")
    return table

def _rich_subcommand_tables(parser):
    action = _subparser_action(parser)
    if action is None:
        return []

    rows_by_name = {}
    for subaction in action._get_subactions():
        command_name = str(subaction.metavar)
        rows_by_name[command_name] = (command_name, subaction.help or "")

    tables = []
    grouped_names = set()
    for title, names in getattr(parser, "_help_command_groups", []):
        rows = []
        for name in names:
            matched_name = None
            if name in rows_by_name:
                matched_name = name
            else:
                for row_name in rows_by_name:
                    if row_name.startswith(f"{name} "):
                        matched_name = row_name
                        break
            if matched_name is not None:
                rows.append(rows_by_name[matched_name])
                grouped_names.add(matched_name)
        table = _rich_command_table(title, rows)
        if table is not None:
            tables.append(table)

    remaining_rows = [rows_by_name[name] for name in rows_by_name if name not in grouped_names]
    remaining_table = _rich_command_table("Commands", remaining_rows)
    if remaining_table is not None:
        tables.append(remaining_table)
    return tables

def _rich_examples_panel(parser):
    examples = getattr(parser, "_help_examples", None) or []
    if not examples:
        return None

    body = RichText()
    for index, example in enumerate(examples):
        if index:
            body.append("\n")
        body.append("$ ", style="drakkar.heading")
        body.append(example, style="drakkar.help")
    return RichPanel(
        body,
        title=RichText(" Examples ", style="drakkar.heading"),
        border_style="drakkar.rule",
        padding=(1, 2),
    )

def _rich_action_tables(parser):
    custom_sections = getattr(parser, "_help_sections", None) or []
    if not custom_sections:
        return [
            table
            for action_group in parser._action_groups
            for table in [_rich_action_table(parser, action_group)]
            if table is not None
        ]

    tables = []
    assigned = set()
    for title, keys in custom_sections:
        actions = []
        for key in keys:
            action = _find_action_by_key(parser, key)
            if (
                action is None
                or action in assigned
                or action.help is argparse.SUPPRESS
                or isinstance(action, argparse._SubParsersAction)
            ):
                continue
            actions.append(action)
            assigned.add(action)
        if not actions:
            continue
        table = _rich_table(title)
        for action in actions:
            option = ", ".join(action.option_strings) if action.option_strings else action.dest
            table.add_row(option, _action_value(action), _expanded_action_help(parser, action))
        tables.append(table)

    remaining = [
        action
        for action in parser._actions
        if action not in assigned
        and action.help is not argparse.SUPPRESS
        and not isinstance(action, argparse._SubParsersAction)
    ]
    if remaining:
        table = _rich_table("General")
        for action in remaining:
            option = ", ".join(action.option_strings) if action.option_strings else action.dest
            table.add_row(option, _action_value(action), _expanded_action_help(parser, action))
        tables.append(table)
    return tables

def _rich_help_renderable(parser):
    if not _rich_available():
        return None

    header = RichText()
    if parser.description:
        header.append(parser.description, style="drakkar.text")
        header.append("\n\n")
    category = getattr(parser, "_help_category", None)
    if category:
        header.append("Category\n", style="drakkar.heading")
        header.append(str(category), style="drakkar.help")
        header.append("\n\n")
    header.append("Usage\n", style="drakkar.heading")
    header.append(_clean_usage(parser), style="drakkar.help")

    renderables = [
        RichPanel(
            header,
            title=RichText(f" {parser.prog} ", style="drakkar.heading"),
            border_style="drakkar.rule",
            padding=(1, 2),
        )
    ]

    renderables.extend(_rich_subcommand_tables(parser))

    examples_panel = _rich_examples_panel(parser)
    if examples_panel is not None:
        renderables.append(examples_panel)

    renderables.extend(_rich_action_tables(parser))

    return RichGroup(*renderables)

class RichArgumentParser(argparse.ArgumentParser):
    """Render argparse output through the shared Rich console."""

    def print_help(self, file=None):
        renderable = _rich_help_renderable(self)
        target = get_console(file)
        if target is not None and renderable is not None:
            display_banner_sequence(get_drakkar_banner_renderables(include_intro=False), file=file, delay_after=True)
            target.print(renderable)
            return
        output_file = file if file is not None else sys.stdout
        display_banner_sequence(get_drakkar_banner_renderables(include_intro=False), file=output_file, delay_after=True)
        super().print_help(output_file)
        if not _rich_available():
            output_file.write(
                "\nNote: Rich-styled help is unavailable because the 'rich' package is not installed in this environment.\n"
                "Reinstall Drakkar with dependencies or run 'python -m pip install rich'.\n"
            )

    def _print_message(self, message, file=None):
        if message:
            print(message, end="", file=file, style="drakkar.help")

    def error(self, message):
        self.print_usage(sys.stderr)
        print(f"ERROR: {message}", file=sys.stderr)
        self.exit(2)
