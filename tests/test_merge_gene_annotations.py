from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "merge_gene_annotations.py"


def load_merge_module():
    spec = importlib.util.spec_from_file_location("merge_gene_annotations", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class MergeGeneAnnotationTests(unittest.TestCase):
    def test_default_identity_threshold_is_50(self) -> None:
        module = load_merge_module()

        self.assertEqual(module.DEFAULT_IDENTITY_THRESHOLD, 50.0)

    def test_vfdb_parser_filters_by_evalue_and_identity(self) -> None:
        # Identity values are in percentage (0-100), matching mmseqs2 pident format.
        # The vfdb rule uses --format-output pident to ensure this range.
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vf_hits = tmp / "vfdb.txt"
            vf_map = tmp / "vfdb.tsv"

            vf_hits.write_text(
                "\n".join(
                    [
                        "gene1\tentry_low_identity\t95\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200",
                        "gene1\tentry_keep\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t180",
                        "gene2\tentry_bad_evalue\t99\t100\t0\t0\t1\t100\t1\t100\t1e-5\t120",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            vf_map.write_text(
                "\n".join(
                    [
                        "entry\tvf\tvf_type",
                        "entry_low_identity\tlow\tlow_type",
                        "entry_keep\tkeep\tkeep_type",
                        "entry_bad_evalue\tbad\tbad_type",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            parsed = module.parse_vfdb(vf_hits, vf_map, 1e-10, 98)

        self.assertEqual(parsed.to_dict("records"), [{"gene": "gene1", "vf": "keep", "vf_type": "keep_type"}])

    def test_vfdb_parser_rejects_fractional_identity(self) -> None:
        # mmseqs2 default output uses fident (0-1) but the vfdb rule explicitly
        # requests pident (0-100). This test confirms that fractional values (0-1)
        # are all filtered out by the percentage threshold, catching any regression
        # where the --format-output flag is accidentally removed.
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vf_hits = tmp / "vfdb.txt"
            vf_map = tmp / "vfdb.tsv"

            vf_hits.write_text(
                "\n".join(
                    [
                        "gene1\tentry_a\t0.99\t100\t1\t0\t1\t100\t1\t100\t1e-50\t200",
                        "gene2\tentry_b\t0.75\t100\t25\t0\t1\t100\t1\t100\t1e-20\t150",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            vf_map.write_text(
                "\n".join(["entry\tvf\tvf_type", "entry_a\tA\ttype_a", "entry_b\tB\ttype_b"])
                + "\n",
                encoding="utf-8",
            )

            parsed = module.parse_vfdb(vf_hits, vf_map, 1e-10, 50.0)

        self.assertEqual(len(parsed), 0, "Fractional identity values must not pass the percentage threshold")


    def test_uniprot_accession_from_target(self) -> None:
        module = load_merge_module()

        self.assertEqual(module.uniprot_accession_from_target("AF-P12345-F1-model_v4.cif.gz"), "P12345")
        self.assertEqual(module.uniprot_accession_from_target("AF-A0A0B1-F1-model_v4"), "A0A0B1")
        # Non-AlphaFold target ids are returned unchanged so they can still join.
        self.assertEqual(module.uniprot_accession_from_target("1abc_A"), "1abc_A")

    def test_foldseek_parser_picks_best_hit_and_maps_function(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            m8 = tmp / "foldseek.m8"
            mapping = tmp / "map.tsv"

            m8.write_text(
                "\n".join(
                    [
                        # gene1 has two hits; the lower e-value (P_keep) must win.
                        "gene1\tAF-P_worse-F1-model_v4.cif.gz\t0.4\t100\t10\t0\t1\t100\t1\t100\t1e-12\t90",
                        "gene1\tAF-P_keep-F1-model_v4.cif.gz\t0.6\t100\t5\t0\t1\t100\t1\t100\t1e-30\t150",
                        # gene2's only hit is above the e-value cutoff -> dropped.
                        "gene2\tAF-P_bad-F1-model_v4.cif.gz\t0.3\t100\t40\t0\t1\t100\t1\t100\t1e-5\t40",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            mapping.write_text(
                "\n".join(
                    [
                        "accession\tkegg\tec\tpfam",
                        "P_keep\tK00010\t1.1.1.1\tPF00010",
                        "P_worse\tK99999\t9.9.9.9\tPF99999",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            parsed = module.parse_foldseek(m8, mapping, 1e-10)

        self.assertEqual(
            parsed.to_dict("records"),
            [{"gene": "gene1", "kegg": "K00010", "ec": "1.1.1.1", "pfam": "PF00010"}],
        )

    def test_structure_fills_only_empty_columns_and_sets_evidence(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            gff = tmp / "genes.gff"
            kegg = tmp / "kegg.tsv"
            m8 = tmp / "foldseek.m8"
            mapping = tmp / "map.tsv"
            out = tmp / "out.tsv"

            gff.write_text(
                "##gff-version 3\n"
                "c1\tProdigal\tCDS\t1\t90\t.\t+\t0\tID=1_1;partial=00\n"
                "c1\tProdigal\tCDS\t100\t200\t.\t+\t0\tID=1_2;partial=00\n"
                "c1\tProdigal\tCDS\t300\t450\t.\t+\t0\tID=1_3;partial=00\n",
                encoding="utf-8",
            )
            # c1_1 gets a sequence KO; structure must not overwrite it.
            kegg.write_text(
                "#columns\n"
                "K00001               -          c1_1        -           1e-30  100.0   0.0   "
                "1e-30  100.0   0.0   1.0   1   0   0   1   1   1   1 desc\n",
                encoding="utf-8",
            )
            # Structure hits both c1_1 (already annotated) and c1_2 (orphan).
            m8.write_text(
                "c1_1\tAF-P_seq-F1-model_v4.cif.gz\t0.5\t100\t10\t0\t1\t100\t1\t100\t1e-20\t120\n"
                "c1_2\tAF-P_str-F1-model_v4.cif.gz\t0.5\t100\t10\t0\t1\t100\t1\t100\t1e-25\t140\n",
                encoding="utf-8",
            )
            mapping.write_text(
                "accession\tkegg\tec\tpfam\n"
                "P_seq\tK00002\t2.2.2.2\tPF00002\n"
                "P_str\tK00003\t3.3.3.3\tPF00003\n",
                encoding="utf-8",
            )

            module.merge_annotations(
                str(gff), str(kegg), "", "", "", "", "", "", "", "", "", "",
                str(out), foldseek_file=str(m8), foldseekdb_file=str(mapping),
            )

            import csv

            rows = {row["gene"]: row for row in csv.DictReader(out.open(), delimiter="\t")}

        # Sequence KO preserved; structure did not overwrite it.
        self.assertEqual(rows["c1_1"]["kegg"], "K00001")
        self.assertEqual(rows["c1_1"]["evidence"], "sequence")
        # Orphan filled by structure.
        self.assertEqual(rows["c1_2"]["kegg"], "K00003")
        self.assertEqual(rows["c1_2"]["ec"], "3.3.3.3")
        self.assertEqual(rows["c1_2"]["evidence"], "structure")
        # No annotation at all.
        self.assertEqual(rows["c1_3"]["evidence"], "")


if __name__ == "__main__":
    unittest.main()
