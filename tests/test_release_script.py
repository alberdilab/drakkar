from __future__ import annotations

from datetime import date
import contextlib
import importlib.util
import io
from pathlib import Path
import re
import sys
import unittest


def _load_release_module():
    module_path = Path("scripts/release.py")
    spec = importlib.util.spec_from_file_location("drakkar_release_script", module_path)
    if spec is None or spec.loader is None:
        raise AssertionError("Could not load scripts/release.py for testing.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class ReleaseScriptTests(unittest.TestCase):
    def test_package_declares_actual_minimum_python(self) -> None:
        pyproject = Path("pyproject.toml").read_text(encoding="utf-8")
        match = re.search(r'^requires-python = "([^"]+)"$', pyproject, re.MULTILINE)
        self.assertIsNotNone(match)
        self.assertEqual(match.group(1), ">=3.10")

    def test_release_workflow_gates_artifacts_on_all_supported_pythons(self) -> None:
        workflow = Path(".github/workflows/release.yml").read_text(encoding="utf-8")
        self.assertIn('python-version: ["3.10", "3.11", "3.12"]', workflow)
        self.assertRegex(workflow, r"(?m)^  build:\n(?:.*\n)*?    needs: test$")
        self.assertIn('test "v${package_version}" = "${GITHUB_REF_NAME}"', workflow)

    def test_release_next_steps_stage_new_files_before_committing(self) -> None:
        release = _load_release_module()
        plan = release.ReleasePlan(
            version="2.0.0",
            release_date="2026-08-24",
            run_tests=True,
            run_build=True,
            run_twine_check=True,
            dry_run=False,
        )
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            release.print_release_plan(plan, [])

        # The actionable next steps are printed by main after the plan. Keep a
        # source-level guard so new files cannot silently be missed again.
        source = Path("scripts/release.py").read_text(encoding="utf-8")
        self.assertIn("git add -A", source)
        self.assertNotIn("git commit -am", source)
        self.assertIn("Wait for the Python 3.10-3.12 test workflow", source)

    def test_update_pyproject_version_rewrites_single_version_line(self) -> None:
        release = _load_release_module()
        content = '[project]\nname = "drakkar"\nversion = "1.0.1"\n'
        updated = release.update_pyproject_version(content, "1.0.2")
        self.assertIn('version = "1.0.2"', updated)
        self.assertNotIn('version = "1.0.1"', updated)

    def test_parse_args_defaults_release_date_to_today(self) -> None:
        release = _load_release_module()
        parsed = release.parse_args(["1.0.2"])
        self.assertEqual(parsed.version, "1.0.2")
        self.assertEqual(parsed.release_date, date.today().isoformat())
        self.assertFalse(parsed.skip_tests)

    def test_update_package_version_rewrites_single_dunder_version_line(self) -> None:
        release = _load_release_module()
        content = '__version__ = "1.0.1"\n'
        updated = release.update_package_version(content, "1.0.2")
        self.assertIn('__version__ = "1.0.2"', updated)
        self.assertNotIn('__version__ = "1.0.1"', updated)

    def test_release_changelog_moves_unreleased_section_into_new_version(self) -> None:
        release = _load_release_module()
        content = """# Changelog

## [Unreleased]

### Added

- New release item

### Fixed

- Important bug fix

## [1.0.1] - 2026-04-21

### Changed

- Previous release item
"""
        updated = release.release_changelog(content, "1.0.2", "2026-04-22")
        self.assertIn("## [Unreleased]", updated)
        self.assertIn("- No unreleased changes yet.", updated)
        self.assertIn("## [1.0.2] - 2026-04-22", updated)
        self.assertIn("- New release item", updated)
        self.assertIn("- Important bug fix", updated)
        self.assertIn("## [1.0.1] - 2026-04-21", updated)
        self.assertIn("## [1.0.2] - 2026-04-22\n\n### Added", updated)
        self.assertIn("- Important bug fix\n\n## [1.0.1] - 2026-04-21", updated)

    def test_release_changelog_rejects_placeholder_only_unreleased_section(self) -> None:
        release = _load_release_module()
        content = """# Changelog

## [Unreleased]

### Added

- No unreleased changes yet.
"""
        with self.assertRaises(release.ReleaseError):
            release.release_changelog(content, "1.0.2", "2026-04-22")

    def test_release_changelog_normalizes_blank_line_between_versions(self) -> None:
        release = _load_release_module()
        content = """# Changelog

## [Unreleased]

### Changed

- New formatting safeguard
## [1.0.1] - 2026-04-21

### Changed

- Previous release item
"""
        updated = release.release_changelog(content, "1.0.2", "2026-04-22")
        self.assertIn("## [1.0.2] - 2026-04-22\n\n### Changed", updated)
        self.assertIn("- New formatting safeguard\n\n## [1.0.1] - 2026-04-21", updated)


if __name__ == "__main__":
    unittest.main()
