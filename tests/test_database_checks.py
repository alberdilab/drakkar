from __future__ import annotations

import contextlib
import io
import tempfile
import types
import unittest
from pathlib import Path

import yaml

from drakkar.database_checks import (
    MANAGED_REQUIRED_ARTIFACTS,
    annotating_requirements,
    check_database_artifacts,
    check_database_provenance,
    collect_database_provenance,
    compare_database_provenance,
    manifest_provenance,
    missing_artifacts,
    module_requirements,
    previous_database_provenance,
    reinstall_command,
    stale_outputs_for,
)


def install_kegg_release(release_dir, *, skip=()):
    release_dir.mkdir(parents=True, exist_ok=True)
    for suffix in MANAGED_REQUIRED_ARTIFACTS["kegg"]:
        if suffix in skip:
            continue
        (release_dir / f"kofams{suffix}").write_text("content", encoding="utf-8")
    return release_dir


def write_run_metadata(output_dir, run_id, databases, legacy=False):
    """Record a run's database provenance.

    ``legacy`` puts the metadata in the output root, where runs predating the
    logging directory wrote it.
    """
    payload = {
        "run_id": run_id,
        "command": "annotating",
        "databases": databases,
    }
    directory = Path(output_dir) if legacy else Path(output_dir) / "logging"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"drakkar_{run_id}.yaml"
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return path


class RequirementSelectionTests(unittest.TestCase):
    def test_annotating_requirements_cover_requested_components_only(self) -> None:
        config = {
            "KEGG_DB": "/db/kofams/2026-02-01",
            "PFAM_DB": "/db/pfam/Pfam37.4",
            "GTDB_DB": "/db/gtdb/release232",
        }
        requirements = annotating_requirements("kegg,taxonomy", config=config)
        self.assertEqual(
            [requirement.config_key for requirement in requirements],
            ["GTDB_DB", "KEGG_DB"],
        )

    def test_annotating_requirements_resolve_requested_gtdb_release(self) -> None:
        config = {"GTDB_DB": "/db/gtdb/default", "GTDB_DB_226": "/db/gtdb/release226"}
        requirements = annotating_requirements("taxonomy", gtdb_version="226", config=config)
        self.assertEqual(requirements[0].configured, "/db/gtdb/release226")

    def test_annotating_requirements_skip_unconfigured_databases(self) -> None:
        self.assertEqual(annotating_requirements("kegg", config={"KEGG_DB": ""}), [])

    def test_module_requirements_include_singlem_only_with_fraction(self) -> None:
        config = {"SINGLEM_DB": "/db/singlem", "CHECKM2_DB": "/db/checkm2.dmnd"}
        args = types.SimpleNamespace(fraction=False)
        self.assertEqual(
            [requirement.config_key for requirement in module_requirements("preprocessing", args, config)],
            [],
        )
        args = types.SimpleNamespace(fraction=True)
        self.assertEqual(
            [requirement.config_key for requirement in module_requirements("preprocessing", args, config)],
            ["SINGLEM_DB"],
        )


class ArtifactCheckTests(unittest.TestCase):
    def test_complete_release_reports_no_missing_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            release = install_kegg_release(Path(tmpdir) / "kofams" / "2026-02-01")
            requirement = annotating_requirements("kegg", config={"KEGG_DB": str(release)})[0]
            self.assertEqual(missing_artifacts(requirement), [])

    def test_missing_ko_list_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            release = install_kegg_release(Path(tmpdir) / "kofams" / "2026-02-01", skip=("_ko_list.tsv",))
            requirement = annotating_requirements("kegg", config={"KEGG_DB": str(release)})[0]
            self.assertEqual(missing_artifacts(requirement), [str(release / "kofams_ko_list.tsv")])

    def test_empty_artifact_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            release = install_kegg_release(Path(tmpdir) / "kofams" / "2026-02-01")
            (release / "kofams.json").write_text("", encoding="utf-8")
            requirement = annotating_requirements("kegg", config={"KEGG_DB": str(release)})[0]
            self.assertEqual(missing_artifacts(requirement), [str(release / "kofams.json")])

    def test_absent_release_directory_reports_the_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            release = Path(tmpdir) / "kofams" / "2026-04-01"
            requirement = annotating_requirements("kegg", config={"KEGG_DB": str(release)})[0]
            self.assertEqual(missing_artifacts(requirement), [str(release)])

    def test_unmanaged_database_only_needs_a_populated_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            models = Path(tmpdir) / "defensefinder"
            requirement = annotating_requirements("defense", config={"DEFENSEFINDER_DB": str(models)})[0]
            self.assertEqual(missing_artifacts(requirement), [str(models)])
            models.mkdir()
            self.assertEqual(missing_artifacts(requirement), [str(models)])
            (models / "models").write_text("content", encoding="utf-8")
            self.assertEqual(missing_artifacts(requirement), [])

    def test_check_database_artifacts_blocks_and_names_the_reinstall_command(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            release = install_kegg_release(Path(tmpdir) / "kofams" / "2026-02-01", skip=("_ko_list.tsv",))
            requirements = annotating_requirements("kegg", config={"KEGG_DB": str(release)})
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                allowed = check_database_artifacts(requirements)
            self.assertFalse(allowed)
            self.assertIn("kofams_ko_list.tsv", stream.getvalue())
            self.assertIn("drakkar database kegg", stream.getvalue())

    def test_check_database_artifacts_can_be_skipped(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            release = Path(tmpdir) / "kofams" / "2026-02-01"
            requirements = annotating_requirements("kegg", config={"KEGG_DB": str(release)})
            self.assertTrue(check_database_artifacts(requirements, skip=True))

    def test_reinstall_command_uses_the_release_directory(self) -> None:
        requirement = annotating_requirements("kegg", config={"KEGG_DB": "/db/kofams/2026-02-01"})[0]
        self.assertEqual(
            reinstall_command(requirement),
            "drakkar database kegg --directory /db/kofams --version 2026-02-01",
        )

    def test_reinstall_command_is_absent_for_unmanaged_databases(self) -> None:
        requirement = annotating_requirements("mobile", config={"GENOMAD_DB": "/db/genomad"})[0]
        self.assertIsNone(reinstall_command(requirement))


class ProvenanceTests(unittest.TestCase):
    def test_collect_provenance_reads_the_installed_version_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            release = install_kegg_release(Path(tmpdir) / "kofams" / "2026-02-01")
            (release / "database_versions.yaml").write_text(
                yaml.safe_dump({"requested_version": "2026-02-01", "source_version": "kofam archive 2026-02-01"}),
                encoding="utf-8",
            )
            requirements = annotating_requirements("kegg", config={"KEGG_DB": str(release)})
            provenance = collect_database_provenance(requirements)
            self.assertEqual(provenance["KEGG_DB"]["requested_version"], "2026-02-01")
            self.assertEqual(provenance["KEGG_DB"]["release"], "2026-02-01")

    def test_collect_provenance_normalizes_a_configured_target_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            release = install_kegg_release(Path(tmpdir) / "kofams" / "2026-02-01")
            requirements = annotating_requirements("kegg", config={"KEGG_DB": str(release / "kofams")})
            provenance = collect_database_provenance(requirements)
            self.assertEqual(provenance["KEGG_DB"]["configured"], str(release))

    def test_compare_ignores_databases_the_earlier_run_did_not_record(self) -> None:
        current = {"KEGG_DB": {"configured": "/db/kofams/2026-04-01", "release": "2026-04-01"}}
        previous = {"PFAM_DB": {"configured": "/db/pfam/Pfam37.4", "release": "Pfam37.4"}}
        self.assertEqual(compare_database_provenance(current, previous), [])

    def test_compare_detects_a_changed_release(self) -> None:
        current = {"KEGG_DB": {"configured": "/db/kofams/2026-04-01", "release": "2026-04-01"}}
        previous = {"KEGG_DB": {"configured": "/db/kofams/2026-02-01", "release": "2026-02-01"}}
        changes = compare_database_provenance(current, previous)
        self.assertEqual([change[0] for change in changes], ["KEGG_DB"])

    def test_reinstalling_the_same_release_is_not_a_change(self) -> None:
        current = {"KEGG_DB": {"configured": "/db/kofams/2026-02-01", "release": "2026-02-01"}}
        previous = {
            "KEGG_DB": {
                "configured": "/db/kofams/2026-02-01",
                "release": "2026-02-01",
                "source_version": "kofam archive 2026-02-01",
            }
        }
        self.assertEqual(compare_database_provenance(current, previous), [])

    def test_previous_provenance_prefers_the_most_recent_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            write_run_metadata(tmpdir, "20260101-000000", {"KEGG_DB": {"configured": "/db/a", "release": "a"}})
            write_run_metadata(tmpdir, "20260201-000000", {"KEGG_DB": {"configured": "/db/b", "release": "b"}})
            source, previous = previous_database_provenance(tmpdir)
            self.assertIn("20260201-000000", source)
            self.assertEqual(previous["KEGG_DB"]["release"], "b")

    def test_previous_provenance_reads_runs_from_the_output_root(self) -> None:
        # A directory written before the logging directory existed still guards
        # against mixing database releases.
        with tempfile.TemporaryDirectory() as tmpdir:
            write_run_metadata(
                tmpdir,
                "20260101-000000",
                {"KEGG_DB": {"configured": "/db/a", "release": "a"}},
                legacy=True,
            )
            source, previous = previous_database_provenance(tmpdir)
            self.assertIn("20260101-000000", source)
            self.assertEqual(previous["KEGG_DB"]["release"], "a")

    def test_previous_provenance_spans_both_layouts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            write_run_metadata(
                tmpdir,
                "20260101-000000",
                {"KEGG_DB": {"configured": "/db/a", "release": "a"}},
                legacy=True,
            )
            write_run_metadata(
                tmpdir, "20260201-000000", {"KEGG_DB": {"configured": "/db/b", "release": "b"}}
            )
            source, previous = previous_database_provenance(tmpdir)
            self.assertIn("20260201-000000", source)
            self.assertEqual(previous["KEGG_DB"]["release"], "b")

    def test_previous_provenance_falls_back_to_the_annotation_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_dir = Path(tmpdir) / "annotating"
            manifest_dir.mkdir()
            (manifest_dir / "annotation_manifest.yaml").write_text(
                yaml.safe_dump(
                    {
                        "databases": {
                            "kegg": {
                                "configured_path": "/db/kofams/2026-02-01",
                                "requested_version": "2026-02-01",
                            }
                        }
                    }
                ),
                encoding="utf-8",
            )
            source, previous = previous_database_provenance(tmpdir)
            self.assertIn("annotation manifest", source)
            self.assertEqual(previous["KEGG_DB"]["requested_version"], "2026-02-01")

    def test_manifest_provenance_is_empty_without_a_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            self.assertEqual(manifest_provenance(tmpdir), (None, {}))


class ProvenanceGateTests(unittest.TestCase):
    def _annotated_output(self, tmpdir):
        output = Path(tmpdir) / "output"
        (output / "annotating" / "kegg").mkdir(parents=True)
        (output / "annotating" / "kegg" / "mag1.tsv").write_text("hits", encoding="utf-8")
        return output

    def test_stale_outputs_are_listed_for_a_changed_database(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self._annotated_output(tmpdir)
            requirement = annotating_requirements("kegg", config={"KEGG_DB": "/db/kofams/2026-04-01"})[0]
            self.assertEqual(
                stale_outputs_for(output, requirement),
                [str(output / "annotating" / "kegg" / "mag1.tsv")],
            )

    def test_changed_database_with_existing_outputs_blocks_the_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self._annotated_output(tmpdir)
            write_run_metadata(
                output,
                "20260101-000000",
                {"KEGG_DB": {"configured": "/db/kofams/2026-02-01", "release": "2026-02-01"}},
            )
            requirements = annotating_requirements("kegg", config={"KEGG_DB": "/db/kofams/2026-04-01"})
            current = collect_database_provenance(requirements)
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                allowed = check_database_provenance(output, requirements, current)
            self.assertFalse(allowed)
            self.assertIn("2026-02-01", stream.getvalue())
            self.assertIn("2026-04-01", stream.getvalue())
            self.assertIn("mag1.tsv", stream.getvalue())

    def test_changed_database_without_existing_outputs_is_informational(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "output"
            output.mkdir()
            write_run_metadata(
                output,
                "20260101-000000",
                {"KEGG_DB": {"configured": "/db/kofams/2026-02-01", "release": "2026-02-01"}},
            )
            requirements = annotating_requirements("kegg", config={"KEGG_DB": "/db/kofams/2026-04-01"})
            current = collect_database_provenance(requirements)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertTrue(check_database_provenance(output, requirements, current))

    def test_allow_database_change_overrides_the_block(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self._annotated_output(tmpdir)
            write_run_metadata(
                output,
                "20260101-000000",
                {"KEGG_DB": {"configured": "/db/kofams/2026-02-01", "release": "2026-02-01"}},
            )
            requirements = annotating_requirements("kegg", config={"KEGG_DB": "/db/kofams/2026-04-01"})
            current = collect_database_provenance(requirements)
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                allowed = check_database_provenance(output, requirements, current, allow_change=True)
            self.assertTrue(allowed)
            self.assertIn("--allow-database-change", stream.getvalue())

    def test_unchanged_databases_pass_silently(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self._annotated_output(tmpdir)
            requirements = annotating_requirements("kegg", config={"KEGG_DB": "/db/kofams/2026-02-01"})
            current = collect_database_provenance(requirements)
            write_run_metadata(output, "20260101-000000", current)
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                self.assertTrue(check_database_provenance(output, requirements, current))
            self.assertEqual(stream.getvalue(), "")

    def test_directory_without_earlier_runs_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self._annotated_output(tmpdir)
            requirements = annotating_requirements("kegg", config={"KEGG_DB": "/db/kofams/2026-04-01"})
            current = collect_database_provenance(requirements)
            self.assertTrue(check_database_provenance(output, requirements, current))


if __name__ == "__main__":
    unittest.main()
