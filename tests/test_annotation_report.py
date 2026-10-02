from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import yaml

import drakkar


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "write_annotation_report.py"


def load_module():
    spec = importlib.util.spec_from_file_location("write_annotation_report", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class AnnotationReportTests(unittest.TestCase):
    def test_manifest_links_qc_tables_tools_and_database_checksums(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            database = tmp / "kegg" / "2026-02-01"
            database.mkdir(parents=True)
            (database / "database_versions.yaml").write_text(
                yaml.safe_dump({
                    "requested_version": "2026-02-01",
                    "source_version": "kofam archive 2026-02-01",
                    "sources": ["https://example.test/kofam"],
                    "files": [{"path": "kofam", "sha256": "abc123", "size_bytes": 42}],
                }),
                encoding="utf-8",
            )
            environment = tmp / "environment.yaml"
            environment.write_text(
                "dependencies:\n  - python=3.12\n  - hmmer=3.4\n",
                encoding="utf-8",
            )
            table = tmp / "gene_annotations.tsv.xz"
            table.write_bytes(b"artifact")
            qc_input = tmp / "MAG_A.qc.json"
            qc_input.write_text(
                json.dumps({
                    "sources": [{
                        "mag": "MAG_A",
                        "level": "gene",
                        "source": "kegg",
                        "reported_records": 5,
                        "retained_records": 3,
                        "rejected_records": 2,
                        "unmapped_records": 0,
                        "unique_entities": 2,
                        "filter_stage": "drakkar",
                    }]
                }),
                encoding="utf-8",
            )
            manifest_output = tmp / "annotation_manifest.yaml"
            qc_output = tmp / "annotation_qc.tsv"

            module.write_annotation_report(
                qc_inputs=[qc_input],
                manifest_output=manifest_output,
                qc_output=qc_output,
                enabled_sources=["kegg"],
                thresholds={"kofam_acceptance": "native_model_cutoffs"},
                databases={"kegg": str(database)},
                tools={"hmmer": "hmmer/3.4"},
                environments={"functional": str(environment)},
                tables={"gene_annotations": str(table)},
            )

            manifest = yaml.safe_load(manifest_output.read_text(encoding="utf-8"))
            with qc_output.open(encoding="utf-8", newline="") as handle:
                qc_rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(manifest["drakkar_version"], drakkar.__version__)
        self.assertEqual(manifest["enabled_sources"], ["kegg"])
        self.assertEqual(manifest["databases"]["kegg"]["files"][0]["sha256"], "abc123")
        self.assertEqual(manifest["configured_tools"]["hmmer"], "hmmer/3.4")
        self.assertEqual(manifest["environments"]["functional"]["dependencies"][1], "hmmer=3.4")
        self.assertEqual(qc_rows[0]["rejected_records"], "2")

    def test_ncbifam_release_and_checksums_reach_manifest_and_qc_sidecar(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            release = tmp / "ncbifam" / "20.0"
            release.mkdir(parents=True)
            files = [
                {
                    "filename": "hmm_PGAP.LIB",
                    "role": "profile_library",
                    "sha256": "library-sha256",
                    "source_url": "https://ftp.ncbi.nlm.nih.gov/hmm/20.0/hmm_PGAP.LIB",
                },
                {
                    "filename": "hmm_PGAP.tsv",
                    "role": "profile_metadata",
                    "sha256": "metadata-sha256",
                    "source_url": "https://ftp.ncbi.nlm.nih.gov/hmm/20.0/hmm_PGAP.tsv",
                },
            ]
            (release / "database_versions.yaml").write_text(
                yaml.safe_dump({
                    "requested_version": "20.0",
                    "source_version": "NCBIfam/PGAP HMM release 20.0",
                    "retrieved_at": "2026-06-25T10:25:00+00:00",
                    "sources": [record["source_url"] for record in files],
                    "files": files,
                }),
                encoding="utf-8",
            )
            qc_input = tmp / "MAG_A.qc.json"
            qc_input.write_text(json.dumps({
                "sources": [{
                    "mag": "MAG_A",
                    "level": "gene",
                    "source": "ncbifam",
                    "reported_records": 2,
                    "retained_records": 2,
                    "rejected_records": 0,
                    "unmapped_records": 0,
                    "unique_entities": 1,
                    "filter_stage": "upstream_native_validated",
                    "database_release": "20.0",
                    "database_source_version": "NCBIfam/PGAP HMM release 20.0",
                    "database_checksums": json.dumps({
                        "hmm_PGAP.LIB": "library-sha256",
                        "hmm_PGAP.tsv": "metadata-sha256",
                    }),
                }]
            }), encoding="utf-8")
            manifest_output = tmp / "annotation_manifest.yaml"
            qc_output = tmp / "annotation_qc.tsv"

            manifest, _ = module.write_annotation_report(
                qc_inputs=[qc_input],
                manifest_output=manifest_output,
                qc_output=qc_output,
                enabled_sources=["ncbifam"],
                thresholds={"ncbifam_acceptance": "native_trusted_cutoffs_no_fallback"},
                databases={"ncbifam": str(release)},
                tools={"hmmer": "hmmer/3.4"},
            )
            with qc_output.open(encoding="utf-8", newline="") as handle:
                qc_row = next(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(manifest["databases"]["ncbifam"]["requested_version"], "20.0")
        self.assertEqual(manifest["databases"]["ncbifam"]["retrieved_at"], "2026-06-25T10:25:00+00:00")
        self.assertEqual(
            {record["filename"]: record["sha256"] for record in manifest["databases"]["ncbifam"]["files"]},
            {"hmm_PGAP.LIB": "library-sha256", "hmm_PGAP.tsv": "metadata-sha256"},
        )
        self.assertEqual(qc_row["database_release"], "20.0")
        self.assertEqual(
            json.loads(qc_row["database_checksums"])["hmm_PGAP.LIB"],
            "library-sha256",
        )


if __name__ == "__main__":
    unittest.main()
