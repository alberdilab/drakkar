from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar import cli as cli_module
from drakkar import cli_parser
from drakkar import environments as environments_module


CURRENT_ENV = """
channels:
  - conda-forge
  - bioconda
dependencies:
  - samtools=1.21
"""

OLD_ENV = """
channels:
  - conda-forge
  - bioconda
dependencies:
  - samtools=1.19
"""


class EnvironmentDirectoryFixture:
    """Builds a conda prefix laid out the way Snakemake deploys environments."""

    def __init__(self, root: str) -> None:
        self.root = Path(root)
        self.env_dir = self.root / "envs"
        self.definitions = self.root / "definitions"
        self.env_dir.mkdir()
        self.definitions.mkdir()
        (self.definitions / "cataloging.yaml").write_text(CURRENT_ENV, encoding="utf-8")

    def deploy(self, env_hash: str, content: str | None, *, files: int = 1) -> Path:
        target = self.env_dir / env_hash
        (target / "conda-meta").mkdir(parents=True)
        for index in range(files):
            (target / f"file{index}.bin").write_bytes(b"x" * 1024)
        if content is not None:
            (self.env_dir / f"{env_hash}.yaml").write_text(content, encoding="utf-8")
        return target

    def scan(self, **kwargs):
        return environments_module.scan_environments(
            self.env_dir, envs_dir=self.definitions, **kwargs
        )


class EnvironmentScanTests(unittest.TestCase):
    def test_matching_definition_is_reported_in_use(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            fixture.deploy("a" * 32, CURRENT_ENV)

            scan = fixture.scan()

            entry = scan["entries"][0]
            self.assertEqual(entry["status"], environments_module.STATUS_IN_USE)
            self.assertEqual(entry["name"], "cataloging.yaml")
            self.assertFalse(entry["removable"])
            self.assertEqual(scan["missing"], [])

    def test_superseded_definition_is_reported_orphan(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            fixture.deploy("b" * 32, OLD_ENV)

            entry = fixture.scan()["entries"][0]

            self.assertEqual(entry["status"], environments_module.STATUS_ORPHAN)
            self.assertIsNone(entry["name"])
            self.assertTrue(entry["removable"])

    def test_formatting_differences_do_not_orphan_an_environment(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            reformatted = CURRENT_ENV.replace("channels:", "channels:  ").replace("\n", "\n\n")
            fixture.deploy("c" * 32, reformatted)

            entry = fixture.scan()["entries"][0]

            self.assertEqual(entry["status"], environments_module.STATUS_IN_USE)

    def test_environment_without_definition_is_incomplete(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            fixture.deploy("d" * 32, None)

            entry = fixture.scan()["entries"][0]

            self.assertEqual(entry["status"], environments_module.STATUS_INCOMPLETE)
            self.assertTrue(entry["removable"])

    def test_non_hash_directories_are_never_removable(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            fixture.deploy("my_personal_env", None)

            entry = fixture.scan()["entries"][0]

            self.assertEqual(entry["status"], environments_module.STATUS_UNKNOWN)
            self.assertFalse(entry["removable"])

    def test_definitions_without_environment_are_leftovers(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            (fixture.env_dir / f"{'e' * 32}.yaml").write_text(OLD_ENV, encoding="utf-8")
            (fixture.env_dir / f"{'e' * 32}.post-deploy.sh").write_text("echo hi", encoding="utf-8")

            scan = fixture.scan()

            self.assertEqual(scan["entries"], [])
            self.assertEqual(len(scan["leftovers"]), 2)

    def test_undeployed_definitions_are_reported_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            (fixture.definitions / "profiling_genomes.yaml").write_text(OLD_ENV, encoding="utf-8")
            fixture.deploy("f" * 32, CURRENT_ENV)

            scan = fixture.scan()

            self.assertEqual(scan["missing"], ["profiling_genomes.yaml"])

    def test_sizes_are_skipped_when_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            fixture.deploy("a" * 32, CURRENT_ENV, files=3)

            self.assertEqual(fixture.scan()["entries"][0]["size"], 3 * 1024)
            self.assertIsNone(fixture.scan(compute_size=False)["entries"][0]["size"])


class EnvironmentPruneTests(unittest.TestCase):
    def run_prune(self, fixture: EnvironmentDirectoryFixture, **kwargs) -> tuple[int, str]:
        buffer = io.StringIO()
        with patch.object(environments_module, "ENVS_DIR", fixture.definitions):
            with contextlib.redirect_stdout(buffer):
                exit_code = environments_module.run_environments_prune(fixture.env_dir, **kwargs)
        return exit_code, buffer.getvalue()

    def test_prune_is_a_dry_run_without_yes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            orphan = fixture.deploy("b" * 32, OLD_ENV)

            exit_code, output = self.run_prune(fixture)

            self.assertEqual(exit_code, 0)
            self.assertIn("Dry run", output)
            self.assertTrue(orphan.exists())

    def test_prune_with_yes_removes_orphans_and_their_definitions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            live = fixture.deploy("a" * 32, CURRENT_ENV)
            orphan = fixture.deploy("b" * 32, OLD_ENV)
            incomplete = fixture.deploy("d" * 32, None)
            personal = fixture.deploy("my_personal_env", None)
            leftover = fixture.env_dir / f"{'e' * 32}.yaml"
            leftover.write_text(OLD_ENV, encoding="utf-8")

            exit_code, _ = self.run_prune(fixture, assume_yes=True)

            self.assertEqual(exit_code, 0)
            self.assertFalse(orphan.exists())
            self.assertFalse((fixture.env_dir / f"{'b' * 32}.yaml").exists())
            self.assertFalse(incomplete.exists())
            self.assertFalse(leftover.exists())
            self.assertTrue(live.exists())
            self.assertTrue((fixture.env_dir / f"{'a' * 32}.yaml").exists())
            self.assertTrue(personal.exists())

    def test_prune_reports_nothing_to_do(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = EnvironmentDirectoryFixture(tmpdir)
            fixture.deploy("a" * 32, CURRENT_ENV)

            exit_code, output = self.run_prune(fixture, assume_yes=True)

            self.assertEqual(exit_code, 0)
            self.assertIn("Nothing to prune", output)

    def test_missing_directory_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                exit_code = environments_module.run_environments_list(Path(tmpdir) / "absent")

            self.assertEqual(exit_code, 1)
            self.assertIn("Environment directory not found", buffer.getvalue())


class EnvironmentCommandTests(unittest.TestCase):
    def parse(self, argv: list[str]):
        return cli_parser.build_parser().parse_args(argv)

    def test_list_flag_dispatches_without_launching_snakemake(self) -> None:
        args = self.parse(["environments", "-e", "/shared/envs", "--list"])
        self.assertTrue(args.list_environments)

        with patch.object(cli_module, "run_environments_list", return_value=0) as listed, \
                patch.object(cli_module, "run_snakemake_environments") as launched, \
                patch.object(cli_module.sys, "argv", ["drakkar", "environments", "-e", "/shared/envs", "--list"]):
            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = cli_module.main()

        self.assertEqual(exit_code, 0)
        launched.assert_not_called()
        listed.assert_called_once_with("/shared/envs", compute_size=True)

    def test_prune_flag_passes_confirmation_and_size_options(self) -> None:
        argv = ["drakkar", "environments", "-e", "/shared/envs", "--prune", "--yes", "--no-size"]

        with patch.object(cli_module, "run_environments_prune", return_value=0) as pruned, \
                patch.object(cli_module, "run_snakemake_environments") as launched, \
                patch.object(cli_module.sys, "argv", argv):
            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = cli_module.main()

        self.assertEqual(exit_code, 0)
        launched.assert_not_called()
        pruned.assert_called_once_with("/shared/envs", assume_yes=True, compute_size=False)

    def test_list_and_prune_are_mutually_exclusive(self) -> None:
        with self.assertRaises(SystemExit):
            with contextlib.redirect_stderr(io.StringIO()):
                self.parse(["environments", "--list", "--prune"])


if __name__ == "__main__":
    unittest.main()
