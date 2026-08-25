from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "merge_cluster_annotations.py"


def load_module():
    spec = importlib.util.spec_from_file_location("merge_cluster_annotations", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class MergeClusterAnnotationTests(unittest.TestCase):
    def test_enabled_sources_exclude_stale_cluster_files_and_preserve_mag(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            dbcan = tmp / "dbcan.tsv"
            stale_antismash = tmp / "antismash.tsv"
            output = tmp / "clusters.tsv"
            qc = tmp / "clusters.qc.json"
            header = (
                "cluster_id\tcontig\tstart\tend\ttype\tgene_count\tsubstrate\t"
                "gene_functions\tpul_id\n"
            )
            dbcan.write_text(
                header + "CGC1\tcontig_1\t10\t200\tPUL\t5\tstarch\tGH13 [2]\tPUL001\n",
                encoding="utf-8",
            )
            stale_antismash.write_text(
                header + "BGC1\tcontig_2\t20\t300\tBGC\t8\tNRPS\tnrps [1]\t\n",
                encoding="utf-8",
            )

            rows = module.merge_cluster_annotations(
                mag="MAG_A",
                output=output,
                sources={"dbcan"},
                dbcan=dbcan,
                antismash=stale_antismash,
                qc_output=qc,
            )
            qc_payload = json.loads(qc.read_text(encoding="utf-8"))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["mag"], "MAG_A")
        self.assertEqual(rows[0]["cluster_id"], "CGC1")
        self.assertEqual(rows[0]["source"], "dbcan")
        self.assertEqual(qc_payload["sources"][0]["unique_entities"], 1)

    def test_defense_system_id_is_not_written_as_contig(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            defense = tmp / "defense.tsv"
            output = tmp / "clusters.tsv"
            defense.write_text(
                "sys_id\treplicon\tsys_beg\tsys_end\ttype\tgenes_count\tsubtype\t"
                "name_of_profiles_in_sys\n"
                "system_1\tcontig_9\t100\t900\tRM\t4\tType_I\tA,B,C,D\n",
                encoding="utf-8",
            )

            rows = module.merge_cluster_annotations(
                mag="MAG_B",
                output=output,
                sources={"defense"},
                defense=defense,
            )

        self.assertEqual(rows[0]["cluster_id"], "system_1")
        self.assertEqual(rows[0]["contig"], "contig_9")
        self.assertEqual(rows[0]["source"], "defensefinder")

    def test_source_schema_changes_fail_loudly(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            broken = tmp / "broken.tsv"
            broken.write_text("contig\tstart\tend\ncontig_1\t1\t2\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing required columns"):
                module.merge_cluster_annotations(
                    mag="MAG_A",
                    output=tmp / "clusters.tsv",
                    sources={"antismash"},
                    antismash=broken,
                )


if __name__ == "__main__":
    unittest.main()
