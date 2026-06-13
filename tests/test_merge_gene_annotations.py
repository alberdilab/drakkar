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


if __name__ == "__main__":
    unittest.main()
