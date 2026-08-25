from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "extract_unannotated_genes.py"


def load_extract_module():
    spec = importlib.util.spec_from_file_location("extract_unannotated_genes", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class ExtractUnannotatedGeneTests(unittest.TestCase):
    def test_dbcan_rows_are_native_accepted_calls_for_structure_gating(self) -> None:
        module = load_extract_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            dbcan = Path(tmpdir) / "dbCAN_hmm_results.tsv"
            dbcan.write_text(
                "HMM Name\tHMM Length\tTarget Name\tTarget Length\ti-Evalue\t"
                "HMM From\tHMM To\tTarget From\tTarget To\tCoverage\tHMM File Name\n"
                "GH5.hmm\t300\tgene1\t500\t1e-20\t1\t150\t20\t170\t0.50\tdbCAN.hmm\n",
                encoding="utf-8",
            )

            annotated = module.annotated_gene_ids([dbcan], 1e-30)

        # The deliberately stricter generic threshold must not override a row
        # that dbCAN has already accepted with its own coverage-aware filters.
        self.assertEqual(annotated, {"gene1"})


if __name__ == "__main__":
    unittest.main()
