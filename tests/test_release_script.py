from __future__ import annotations

from datetime import date
import importlib.util
from pathlib import Path
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

    def test_release_changelog_rejects_placeholder_only_unreleased_section(self) -> None:
        release = _load_release_module()
        content = """# Changelog

## [Unreleased]

### Added

- No unreleased changes yet.
"""
        with self.assertRaises(release.ReleaseError):
            release.release_changelog(content, "1.0.2", "2026-04-22")


if __name__ == "__main__":
    unittest.main()
