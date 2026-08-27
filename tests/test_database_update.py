from __future__ import annotations

import contextlib
import io
import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

from drakkar import cli as cli_module
from drakkar import cli_main
from drakkar import cli_parser
from drakkar import database_update
from drakkar.database_latest import (
    STATUS_AHEAD,
    STATUS_CURRENT,
    STATUS_OUTDATED,
    STATUS_UNCONFIGURED,
    STATUS_UNKNOWN,
    LatestRelease,
)
from drakkar.database_update import (
    apply_database_updates,
    plan_database_updates,
    run_database_update,
)


def release(name, config_key, installed, latest, status, configured, detail=""):
    return LatestRelease(
        name=name,
        label=name.upper(),
        config_key=config_key,
        configured=configured,
        installed=installed,
        latest=latest,
        status=status,
        detail=detail,
    )


PFAM_OUTDATED = release(
    "pfam", "PFAM_DB", "Pfam37.4", "Pfam38.2", STATUS_OUTDATED, "/db/pfam/Pfam37.4"
)
KEGG_OUTDATED = release(
    "kegg", "KEGG_DB", "2026-02-01", "2026-07-02", STATUS_OUTDATED, "/db/kofams/2026-02-01"
)


class PlanTests(unittest.TestCase):
    def plan(self, results):
        with patch.object(database_update, "resolve_all", return_value=results):
            return plan_database_updates(["pfam"], timeout=1)

    def test_outdated_managed_database_becomes_an_update(self) -> None:
        updates, skipped = self.plan([PFAM_OUTDATED])
        self.assertEqual(len(updates), 1)
        self.assertEqual(skipped, [])
        self.assertEqual(updates[0].name, "pfam")
        self.assertEqual(updates[0].latest, "Pfam38.2")

    def test_new_release_is_installed_beside_the_configured_one(self) -> None:
        updates, _ = self.plan([PFAM_OUTDATED])
        self.assertEqual(updates[0].base_directory, Path("/db/pfam"))
        self.assertEqual(updates[0].release_dir, Path("/db/pfam/Pfam38.2"))

    def test_up_to_date_database_is_neither_updated_nor_skipped(self) -> None:
        current = release(
            "pfam", "PFAM_DB", "Pfam38.2", "Pfam38.2", STATUS_CURRENT, "/db/pfam/Pfam38.2"
        )
        updates, skipped = self.plan([current])
        self.assertEqual(updates, [])
        self.assertEqual(skipped, [])

    def test_outdated_database_drakkar_cannot_install_is_skipped(self) -> None:
        gtdb = release("gtdb", "GTDB_DB", "226", "232", STATUS_OUTDATED, "/db/gtdb/release226")
        updates, skipped = self.plan([gtdb])
        self.assertEqual(updates, [])
        self.assertEqual(len(skipped), 1)
        self.assertIn("gtdbtk download-db.sh", skipped[0][1])

    def test_unresolved_database_is_skipped_with_its_reason(self) -> None:
        unknown = release(
            "pfam", "PFAM_DB", None, None, STATUS_UNKNOWN, "/db/pfam/Pfam37.4", "host is down"
        )
        updates, skipped = self.plan([unknown])
        self.assertEqual(updates, [])
        self.assertEqual(skipped[0][1], "host is down")

    def test_unconfigured_database_is_skipped(self) -> None:
        unconfigured = release(
            "pfam", "PFAM_DB", None, "Pfam38.2", STATUS_UNCONFIGURED, None, "PFAM_DB is empty"
        )
        updates, skipped = self.plan([unconfigured])
        self.assertEqual(updates, [])
        self.assertEqual(len(skipped), 1)

    def test_release_ahead_of_the_source_is_skipped_rather_than_downgraded(self) -> None:
        ahead = release(
            "pfam", "PFAM_DB", "Pfam39.0", "Pfam38.2", STATUS_AHEAD, "/db/pfam/Pfam39.0"
        )
        updates, skipped = self.plan([ahead])
        self.assertEqual(updates, [])
        self.assertIn("newer than the source", skipped[0][1])


class ApplyTests(unittest.TestCase):
    def _updates(self):
        with patch.object(database_update, "resolve_all", return_value=[PFAM_OUTDATED, KEGG_OUTDATED]):
            updates, _ = plan_database_updates(["pfam", "kegg"], timeout=1)
        return updates

    def test_every_release_is_installed_in_turn(self) -> None:
        installed = []
        with contextlib.redirect_stdout(io.StringIO()):
            exit_code = apply_database_updates(
                self._updates(), lambda update: installed.append(update.name) or True
            )
        self.assertEqual(exit_code, 0)
        self.assertEqual(installed, ["pfam", "kegg"])

    def test_a_failed_install_does_not_stop_the_others(self) -> None:
        attempted = []

        def install(update):
            attempted.append(update.name)
            return update.name != "pfam"

        with contextlib.redirect_stdout(io.StringIO()) as output:
            exit_code = apply_database_updates(self._updates(), install)

        rendered = output.getvalue()
        self.assertEqual(attempted, ["pfam", "kegg"])
        self.assertEqual(exit_code, 1)
        self.assertIn("1 installed, 1 failed", rendered)
        self.assertIn("drakkar logging -o /db/pfam/Pfam38.2 --failures", rendered)

    def test_leaving_config_untouched_reports_the_paths_to_set(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()) as output:
            apply_database_updates(self._updates(), lambda update: True, set_default=False)
        rendered = output.getvalue()
        self.assertIn("--no-set-default", rendered)
        self.assertIn("PFAM_DB: /db/pfam/Pfam38.2", rendered)

    def test_successful_update_warns_about_outputs_built_with_the_old_release(self) -> None:
        with contextlib.redirect_stdout(io.StringIO()) as output:
            apply_database_updates(self._updates(), lambda update: True)
        self.assertIn("--allow-database-change", output.getvalue())


class RunDatabaseUpdateTests(unittest.TestCase):
    def test_dry_run_reports_the_plan_without_installing_anything(self) -> None:
        installs = []
        with patch.object(database_update, "resolve_all", return_value=[PFAM_OUTDATED]):
            with contextlib.redirect_stdout(io.StringIO()) as output:
                exit_code = run_database_update(
                    ["pfam"], install=lambda update: installs.append(update) or True
                )
        rendered = output.getvalue()
        self.assertEqual(exit_code, 0)
        self.assertEqual(installs, [])
        self.assertIn("Dry run", rendered)
        self.assertIn("Pfam38.2", rendered)

    def test_yes_installs_the_planned_releases(self) -> None:
        installs = []
        with patch.object(database_update, "resolve_all", return_value=[PFAM_OUTDATED]):
            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = run_database_update(
                    ["pfam"],
                    install=lambda update: installs.append(update.name) or True,
                    assume_yes=True,
                )
        self.assertEqual(exit_code, 0)
        self.assertEqual(installs, ["pfam"])

    def test_nothing_to_do_is_reported_without_calling_the_installer(self) -> None:
        current = release(
            "pfam", "PFAM_DB", "Pfam38.2", "Pfam38.2", STATUS_CURRENT, "/db/pfam/Pfam38.2"
        )
        with patch.object(database_update, "resolve_all", return_value=[current]):
            with contextlib.redirect_stdout(io.StringIO()) as output:
                exit_code = run_database_update(["pfam"], install=None, assume_yes=True)
        self.assertEqual(exit_code, 0)
        self.assertIn("0 to install", output.getvalue())

    def test_unknown_database_name_fails_before_any_source_is_queried(self) -> None:
        with patch.object(database_update, "resolve_all") as resolve_all:
            with contextlib.redirect_stdout(io.StringIO()) as output:
                exit_code = run_database_update(["nonsense"])
        self.assertEqual(exit_code, 1)
        resolve_all.assert_not_called()
        self.assertIn("Unknown database", output.getvalue())


class DatabaseUpdateDispatchTests(unittest.TestCase):
    def test_update_is_parsed_as_a_database_subcommand(self) -> None:
        args = cli_parser.build_parser().parse_args(["database", "update", "kegg", "--yes"])
        self.assertEqual(args.database_name, "update")
        self.assertEqual(args.databases, ["kegg"])
        self.assertTrue(args.assume_yes)
        self.assertTrue(args.set_default)

    def test_dry_run_does_not_prompt_about_the_screen_session(self) -> None:
        argv = ["drakkar", "database", "update"]
        with patch.object(cli_module, "run_database_update", return_value=0) as updated, \
                patch.object(cli_module, "run_snakemake_database") as launched, \
                patch.object(cli_main, "check_screen_session") as screen, \
                patch.object(cli_module.sys, "argv", argv):
            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = cli_module.main()

        self.assertEqual(exit_code, 0)
        launched.assert_not_called()
        screen.assert_not_called()
        self.assertEqual(updated.call_args.args[0], [])
        self.assertFalse(updated.call_args.kwargs["assume_yes"])

    def test_installing_warns_about_the_screen_session(self) -> None:
        argv = ["drakkar", "database", "update", "--yes"]
        with patch.object(cli_module, "run_database_update", return_value=0), \
                patch.object(cli_main, "check_screen_session") as screen, \
                patch.object(cli_module.sys, "argv", argv):
            with contextlib.redirect_stdout(io.StringIO()):
                cli_module.main()
        screen.assert_called_once()

    def test_update_does_not_reach_the_single_release_install_path(self) -> None:
        argv = ["drakkar", "database", "update", "--yes"]
        with patch.object(cli_module, "run_database_update", return_value=0), \
                patch.object(cli_module, "run_snakemake_database") as launched, \
                patch.object(cli_main, "check_screen_session"), \
                patch.object(cli_module, "prepare_output_directory") as prepared, \
                patch.object(cli_module.sys, "argv", argv):
            with contextlib.redirect_stdout(io.StringIO()):
                cli_module.main()
        launched.assert_not_called()
        prepared.assert_not_called()

    def test_install_closure_launches_the_workflow_for_each_release(self) -> None:
        """Exercise the real dispatch wiring, not just the entry point."""
        argv = ["drakkar", "database", "update", "pfam", "--yes"]
        with patch.object(database_update, "resolve_all", return_value=[PFAM_OUTDATED]), \
                patch.object(cli_main, "check_screen_session"), \
                patch.object(cli_module, "prepare_output_directory", return_value=True), \
                patch.object(cli_module, "write_launch_metadata", return_value={"run_id": "1"}), \
                patch.object(cli_module, "set_default_database_path", return_value={"PFAM_DB": "/db/pfam/Pfam38.2"}) as defaulted, \
                patch.object(cli_module, "run_snakemake_database") as launched, \
                patch.object(cli_module.sys, "argv", argv):
            with contextlib.redirect_stdout(io.StringIO()) as output:
                exit_code = cli_module.main()

        self.assertEqual(exit_code, 0)
        launched.assert_called_once()
        positional = launched.call_args.args
        self.assertEqual(positional[0], "database")
        self.assertEqual(positional[2], Path("/db/pfam/Pfam38.2"))
        self.assertEqual(positional[5], "pfam")
        self.assertEqual(positional[6], Path("/db/pfam"))
        self.assertEqual(positional[7], "Pfam38.2")
        defaulted.assert_called_once_with("pfam", Path("/db/pfam"), "Pfam38.2")
        self.assertIn("1 installed, 0 failed", output.getvalue())

    def test_a_workflow_failure_is_contained_and_reported(self) -> None:
        argv = ["drakkar", "database", "update", "pfam", "--yes"]
        with patch.object(database_update, "resolve_all", return_value=[PFAM_OUTDATED]), \
                patch.object(cli_main, "check_screen_session"), \
                patch.object(cli_module, "prepare_output_directory", return_value=True), \
                patch.object(cli_module, "write_launch_metadata", return_value={"run_id": "1"}), \
                patch.object(cli_module, "set_default_database_path") as defaulted, \
                patch.object(
                    cli_module,
                    "run_snakemake_database",
                    side_effect=subprocess.CalledProcessError(1, "snakemake"),
                ), \
                patch.object(cli_module.sys, "argv", argv):
            with contextlib.redirect_stdout(io.StringIO()) as output:
                exit_code = cli_module.main()

        self.assertEqual(exit_code, 1)
        # config.yaml must not be repointed at a release that failed to install.
        defaulted.assert_not_called()
        self.assertIn("0 installed, 1 failed", output.getvalue())

    def test_invalid_download_runtime_is_rejected(self) -> None:
        argv = ["drakkar", "database", "update", "--download-runtime", "0"]
        with patch.object(cli_module, "run_database_update") as updated, \
                patch.object(cli_module.sys, "argv", argv):
            with contextlib.redirect_stdout(io.StringIO()):
                cli_module.main()
        updated.assert_not_called()


if __name__ == "__main__":
    unittest.main()
