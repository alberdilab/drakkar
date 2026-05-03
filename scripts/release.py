#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import date
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT_PATH = ROOT / "pyproject.toml"
PACKAGE_INIT_PATH = ROOT / "drakkar" / "__init__.py"
CHANGELOG_PATH = ROOT / "CHANGELOG.md"
UNRELEASED_HEADER = "## [Unreleased]"
UNRELEASED_PLACEHOLDER = "### Added\n\n- No unreleased changes yet.\n"


class ReleaseError(RuntimeError):
    pass


@dataclass(frozen=True)
class ReleasePlan:
    version: str
    release_date: str
    run_tests: bool
    run_build: bool
    run_twine_check: bool
    dry_run: bool


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare a Drakkar release by cutting the CHANGELOG entry from "
            "Unreleased, updating package version metadata, and running release checks."
        )
    )
    parser.add_argument(
        "version",
        help="Release version to write into pyproject.toml, drakkar/__init__.py, and CHANGELOG.md.",
    )
    parser.add_argument(
        "--date",
        dest="release_date",
        default=date.today().isoformat(),
        help="Release date for the changelog entry in YYYY-MM-DD format. Defaults to today.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show the planned changes without modifying files or running commands.",
    )
    parser.add_argument(
        "--skip-tests",
        action="store_true",
        help="Skip running the test suite.",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip building the sdist and wheel artifacts.",
    )
    parser.add_argument(
        "--skip-twine-check",
        action="store_true",
        help="Skip twine validation of built artifacts.",
    )
    return parser.parse_args(argv)


def main() -> int:
    args = parse_args()
    validate_version(args.version)
    validate_release_date(args.release_date)

    plan = ReleasePlan(
        version=args.version,
        release_date=args.release_date,
        run_tests=not args.skip_tests,
        run_build=not args.skip_build,
        run_twine_check=not args.skip_twine_check,
        dry_run=args.dry_run,
    )

    pyproject_before = PYPROJECT_PATH.read_text(encoding="utf-8")
    package_init_before = PACKAGE_INIT_PATH.read_text(encoding="utf-8")
    changelog_before = CHANGELOG_PATH.read_text(encoding="utf-8")

    pyproject_after = update_pyproject_version(pyproject_before, plan.version)
    package_init_after = update_package_version(package_init_before, plan.version)
    changelog_after = release_changelog(changelog_before, plan.version, plan.release_date)

    changed_files: list[str] = []
    if pyproject_before != pyproject_after:
        changed_files.append(PYPROJECT_PATH.relative_to(ROOT).as_posix())
    if package_init_before != package_init_after:
        changed_files.append(PACKAGE_INIT_PATH.relative_to(ROOT).as_posix())
    if changelog_before != changelog_after:
        changed_files.append(CHANGELOG_PATH.relative_to(ROOT).as_posix())

    print_release_plan(plan, changed_files)
    if plan.dry_run:
        return 0

    if pyproject_before != pyproject_after:
        PYPROJECT_PATH.write_text(pyproject_after, encoding="utf-8")
    if package_init_before != package_init_after:
        PACKAGE_INIT_PATH.write_text(package_init_after, encoding="utf-8")
    if changelog_before != changelog_after:
        CHANGELOG_PATH.write_text(changelog_after, encoding="utf-8")

    if plan.run_tests:
        run_command([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"])
    if plan.run_build:
        require_module("build")
        clean_dist()
        run_command([sys.executable, "-m", "build"])
    if plan.run_twine_check:
        require_module("twine")
        run_command([sys.executable, "-m", "twine", "check", "dist/*"], shell=True)

    print()
    print("Release preparation completed.")
    print("Next steps:")
    print("1. Review changes: git diff")
    print(f'2. Commit them: git commit -am "Release v{plan.version}"')
    print(f"3. Tag the release: git tag v{plan.version}")
    print(f"4. Push branch and tag: git push origin main && git push origin v{plan.version}")
    print(f"5. Monitor the GitHub Actions release workflow for tag v{plan.version}")
    print(
        "6. The workflow will build artifacts and create/update the GitHub release. "
        "PyPI publishing runs only when the repository is public and Trusted Publishing is configured."
    )
    return 0


def validate_version(version: str) -> None:
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        raise ReleaseError(
            f"Unsupported version format: {version}. Expected semantic version like 1.2.3."
        )


def validate_release_date(release_date: str) -> None:
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", release_date):
        raise ReleaseError(
            f"Unsupported release date format: {release_date}. Expected YYYY-MM-DD."
        )


def update_pyproject_version(content: str, version: str) -> str:
    updated, count = re.subn(
        r'(?m)^version = "[^"]+"$',
        f'version = "{version}"',
        content,
        count=1,
    )
    if count != 1:
        raise ReleaseError("Could not find a unique version line in pyproject.toml.")
    return updated


def update_package_version(content: str, version: str) -> str:
    updated, count = re.subn(
        r'(?m)^__version__ = "[^"]+"$',
        f'__version__ = "{version}"',
        content,
        count=1,
    )
    if count != 1:
        raise ReleaseError("Could not find a unique __version__ line in drakkar/__init__.py.")
    return updated


def release_changelog(content: str, version: str, release_date: str) -> str:
    if f"## [{version}] -" in content:
        raise ReleaseError(f"CHANGELOG.md already contains an entry for version {version}.")

    unreleased_start = content.find(UNRELEASED_HEADER)
    if unreleased_start == -1:
        raise ReleaseError("CHANGELOG.md does not contain an [Unreleased] section.")

    next_section = content.find("\n## [", unreleased_start + len(UNRELEASED_HEADER))
    if next_section == -1:
        unreleased_block = content[unreleased_start:]
        remainder = ""
    else:
        unreleased_block = content[unreleased_start:next_section]
        remainder = content[next_section:]

    header_line, _, unreleased_body = unreleased_block.partition("\n")
    normalized_body = unreleased_body.strip("\n")
    if not normalized_body.strip():
        raise ReleaseError("The [Unreleased] section is empty. Add release notes before cutting a release.")

    release_body = normalized_body
    if release_body == UNRELEASED_PLACEHOLDER.strip():
        raise ReleaseError(
            "The [Unreleased] section only contains the placeholder entry. "
            "Replace it with real release notes before cutting a release."
        )

    released_section = f"## [{version}] - {release_date}\n\n{release_body}\n"
    new_unreleased = f"{header_line}\n\n{UNRELEASED_PLACEHOLDER.strip()}\n"

    sections = [
        content[:unreleased_start].rstrip("\n"),
        new_unreleased.rstrip("\n"),
        released_section.rstrip("\n"),
    ]
    remainder = remainder.lstrip("\n").rstrip("\n")
    if remainder:
        sections.append(remainder)
    return "\n\n".join(sections) + "\n"


def require_module(module_name: str) -> None:
    result = subprocess.run(
        [sys.executable, "-c", f"import {module_name}"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode == 0:
        return
    raise ReleaseError(
        f"Missing required module '{module_name}'. Install the release tooling with: "
        f"{sys.executable} -m pip install '.[release]'"
    )


def clean_dist() -> None:
    dist_path = ROOT / "dist"
    if dist_path.exists():
        shutil.rmtree(dist_path)


def run_command(command: list[str], *, shell: bool = False) -> None:
    rendered = " ".join(command)
    print(f"Running: {rendered}")
    subprocess.run(
        rendered if shell else command,
        cwd=ROOT,
        shell=shell,
        check=True,
    )


def print_release_plan(plan: ReleasePlan, changed_files: list[str]) -> None:
    print("Drakkar release plan")
    print(f"  version: {plan.version}")
    print(f"  release_date: {plan.release_date}")
    print(f"  dry_run: {plan.dry_run}")
    print(f"  run_tests: {plan.run_tests}")
    print(f"  run_build: {plan.run_build}")
    print(f"  run_twine_check: {plan.run_twine_check}")
    print("  changed_files:")
    if changed_files:
        for path in changed_files:
            print(f"    - {path}")
    else:
        print("    - none")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ReleaseError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
