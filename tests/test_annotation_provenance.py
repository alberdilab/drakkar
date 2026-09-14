from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path

import yaml

from drakkar.database_checks import (
    ANNOTATION_REPORT_SOURCE_NAMES,
    CLUSTER_ANNOTATION_COMPONENTS,
    GENE_ANNOTATION_COMPONENTS,
    annotation_provenance_changes,
    annotation_report_sources,
    check_annotation_provenance,
    read_annotation_manifest,
    stale_annotation_outputs,
)

ROOT = Path(__file__).resolve().parents[1]
ANNOTATION_RULES = ROOT / "drakkar" / "workflow" / "rules" / "annotating_function.smk"


def write_output_dir(
    directory: Path,
    *,
    version: str = "2.5.7",
    sources=("cazy", "genomad", "kegg", "ncbi_amrfinder", "pfam", "signalp", "vfdb"),
    tables=("final/MAG_A_genes.tsv", "final/MAG_A_clusters.tsv", "gene_annotations.tsv.xz"),
) -> Path:
    annotating = directory / "annotating"
    (annotating / "final").mkdir(parents=True, exist_ok=True)
    (annotating / "annotation_manifest.yaml").write_text(
        yaml.safe_dump(
            {
                "schema_version": "drakkar-annotation-manifest-v1",
                "drakkar_version": version,
                "enabled_sources": sorted(sources),
                "thresholds": {},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    for table in tables:
        (annotating / table).touch()
    return directory


class AnnotationProvenanceTests(unittest.TestCase):
    def test_component_lists_match_the_workflow_rules(self) -> None:
        """The CLI-side split must track ENABLED_*_SOURCES in the .smk."""
        rules = ANNOTATION_RULES.read_text(encoding="utf-8")

        def listed(variable: str) -> set[str]:
            block = re.search(rf"{variable} = \[(.*?)\n\]", rules, re.DOTALL)
            assert block is not None, f"{variable} not found in the rules file"
            return set(re.findall(r'\("([a-z]+)",', block.group(1)))

        gene = listed("ENABLED_GENE_SOURCES")
        cluster = listed("ENABLED_CLUSTER_SOURCES")
        # "structure" is work in progress and rejected by the CLI, so it is not
        # a component the provenance check can ever be asked about.
        self.assertEqual(gene - {"structure"}, set(GENE_ANNOTATION_COMPONENTS))
        self.assertEqual(cluster, set(CLUSTER_ANNOTATION_COMPONENTS))

    def test_report_source_names_match_the_workflow_rules(self) -> None:
        rules = ANNOTATION_RULES.read_text(encoding="utf-8")
        block = re.search(r"REPORT_SOURCE_NAMES = \{(.*?)\n\}", rules, re.DOTALL)
        self.assertIsNotNone(block)
        mapping = dict(re.findall(r'"([a-z]+)": "([a-z_]+)",', block.group(1)))
        mapping.pop("structure", None)
        self.assertEqual(mapping, ANNOTATION_REPORT_SOURCE_NAMES)

    def test_annotation_report_sources_uses_manifest_names(self) -> None:
        self.assertEqual(
            annotation_report_sources("kegg,virulence,amr,card,mobile,defense"),
            {"kegg", "vfdb", "ncbi_amrfinder", "card", "genomad", "defensefinder"},
        )
        # Bundle keywords and taxonomy are not annotation sources.
        self.assertEqual(annotation_report_sources("taxonomy,function,genes"), set())

    def test_missing_manifest_is_not_a_change(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            self.assertTrue(check_annotation_provenance(tmpdir, "kegg", "2.6.0"))

    def test_identical_setup_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            write_output_dir(Path(tmpdir))
            self.assertTrue(
                check_annotation_provenance(
                    tmpdir, "kegg,cazy,pfam,virulence,amr,signalp,mobile", "2.5.7"
                )
            )

    def test_added_gene_source_invalidates_only_gene_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            write_output_dir(Path(tmpdir))
            _, manifest = read_annotation_manifest(tmpdir)
            changes, levels = annotation_provenance_changes(
                manifest, "kegg,cazy,pfam,virulence,amr,card,signalp,mobile", "2.5.7"
            )

        self.assertEqual([field for field, _, _ in changes], ["annotation sources"])
        self.assertEqual(levels, {"gene"})

    def test_added_cluster_source_invalidates_only_cluster_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            write_output_dir(Path(tmpdir))
            _, manifest = read_annotation_manifest(tmpdir)
            changes, levels = annotation_provenance_changes(
                manifest, "kegg,cazy,pfam,virulence,amr,signalp,mobile,antismash", "2.5.7"
            )

        self.assertEqual(levels, {"cluster"})

    def test_defense_invalidates_both_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            write_output_dir(Path(tmpdir))
            _, manifest = read_annotation_manifest(tmpdir)
            _, levels = annotation_provenance_changes(
                manifest, "kegg,cazy,pfam,virulence,amr,signalp,mobile,defense", "2.5.7"
            )

        self.assertEqual(levels, {"gene", "cluster"})

    def test_version_change_alone_invalidates_both_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            write_output_dir(Path(tmpdir))
            _, manifest = read_annotation_manifest(tmpdir)
            changes, levels = annotation_provenance_changes(
                manifest, "kegg,cazy,pfam,virulence,amr,signalp,mobile", "2.6.0"
            )

        self.assertEqual([field for field, _, _ in changes], ["drakkar version"])
        self.assertEqual(levels, {"gene", "cluster"})

    def test_the_2_5_to_2_6_amr_upgrade_is_blocked(self) -> None:
        """The upgrade that motivated this check: amr changed database and card appeared."""
        with tempfile.TemporaryDirectory() as tmpdir:
            write_output_dir(Path(tmpdir))
            blocked = not check_annotation_provenance(
                tmpdir, "kegg,cazy,pfam,virulence,amr,card,signalp,mobile", "2.6.0"
            )

        self.assertTrue(blocked)

    def test_change_without_existing_tables_is_informational(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            write_output_dir(Path(tmpdir), tables=())
            self.assertTrue(
                check_annotation_provenance(
                    tmpdir, "kegg,cazy,pfam,virulence,amr,card,signalp,mobile", "2.6.0"
                )
            )

    def test_allow_annotation_change_continues(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            write_output_dir(Path(tmpdir))
            self.assertTrue(
                check_annotation_provenance(
                    tmpdir,
                    "kegg,cazy,pfam,virulence,amr,card,signalp,mobile",
                    "2.6.0",
                    allow_change=True,
                )
            )

    def test_unreadable_manifest_does_not_block(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            annotating = Path(tmpdir) / "annotating"
            annotating.mkdir(parents=True)
            (annotating / "annotation_manifest.yaml").write_text(": not yaml :", encoding="utf-8")
            self.assertTrue(check_annotation_provenance(tmpdir, "kegg", "2.6.0"))

    def test_stale_outputs_lists_only_the_affected_level(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            write_output_dir(Path(tmpdir))
            gene = stale_annotation_outputs(tmpdir, {"gene"})
            cluster = stale_annotation_outputs(tmpdir, {"cluster"})

        self.assertEqual(
            {Path(path).name for path in gene},
            {"MAG_A_genes.tsv", "gene_annotations.tsv.xz"},
        )
        self.assertEqual({Path(path).name for path in cluster}, {"MAG_A_clusters.tsv"})


if __name__ == "__main__":
    unittest.main()
