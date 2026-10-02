from __future__ import annotations

import csv
import importlib.util
import json
import lzma
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "drakkar" / "workflow" / "scripts" / "project_gifter_input.py"


def load_module():
    spec = importlib.util.spec_from_file_location("project_gifter_input", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class GifterInputProjectionTests(unittest.TestCase):
    def test_projection_is_gene_resolved_generic_and_deterministic(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source = tmp / "gene_annotations.tsv.xz"
            output = tmp / "gifter_input.tsv.xz"
            columns = ["mag", "gene", "source", "annotation_id", "details"]
            rows = [
                ["MAG_B", "g1", "prodigal", "CDS", "{}"],
                ["MAG_B", "g1", "kegg", "K00002", json.dumps({"ec": "2.2.2.2;1.1.1.1"})],
                ["MAG_B", "g1", "pfam", "PF00001", json.dumps({
                    "ec_associations": [
                        {"ec": "3.3.3.3", "confidence": 0.9},
                        {"ec": "1.1.1.1", "confidence": 0.8},
                    ]
                })],
                ["MAG_B", "g1", "cazy", "GH5", "{}"],
                # Repeated native domains do not create duplicate marker rows.
                ["MAG_B", "g1", "cazy", "GH5", "{}"],
                ["MAG_B", "g2", "ncbifam", "NF040708.3", json.dumps({
                    "profile_metadata": {"ec_numbers": "4.1.1.111,5.4.3.2"}
                })],
                ["MAG_B", "g2", "tigrfam", "TIGR02053", "{}"],
                # Unrelated Drakkar evidence is deliberately absent from the projection.
                ["MAG_B", "g2", "vfdb", "VFG0001", "{}"],
            ]
            with lzma.open(source, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(columns)
                writer.writerows(rows)

            count = module.project_gifter_input(source, output)
            with lzma.open(output, "rt", encoding="utf-8", newline="") as handle:
                projected = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(count, 10)
        self.assertEqual(list(projected[0]), module.OUTPUT_COLUMNS)
        self.assertEqual(
            [(row["genome_id"], row["gene_id"], row["namespace"], row["accession"]) for row in projected],
            [
                ("MAG_B", "g1", "CAZY", "GH5"),
                ("MAG_B", "g1", "EC", "1.1.1.1"),
                ("MAG_B", "g1", "EC", "2.2.2.2"),
                ("MAG_B", "g1", "EC", "3.3.3.3"),
                ("MAG_B", "g1", "KO", "K00002"),
                ("MAG_B", "g1", "PFAM", "PF00001"),
                ("MAG_B", "g2", "EC", "4.1.1.111"),
                ("MAG_B", "g2", "EC", "5.4.3.2"),
                ("MAG_B", "g2", "NCBIFAM", "NF040708.3"),
                ("MAG_B", "g2", "TIGRFAM", "TIGR02053"),
            ],
        )

    def test_projection_rejects_invalid_details_instead_of_silently_dropping_evidence(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source = tmp / "genes.tsv"
            output = tmp / "gifter.tsv"
            source.write_text(
                "mag\tgene\tsource\tannotation_id\tdetails\n"
                "MAG_A\tg1\tkegg\tK00001\t{broken}\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "Invalid annotation details JSON"):
                module.project_gifter_input(source, output)


if __name__ == "__main__":
    unittest.main()
