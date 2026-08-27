import argparse
import csv
import importlib.util
import json
import lzma
import sys
import tempfile
import unittest
from pathlib import Path

import yaml


SCRIPTS = Path(__file__).resolve().parents[1] / "drakkar" / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location("aggregate_amr", SCRIPTS / "aggregate_amr.py")
aggregate_amr = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(aggregate_amr)


class AmrAggregateTests(unittest.TestCase):
    def test_zero_hit_assemblies_are_kept_in_summary_and_manifest(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_manifest = root / "data" / "amr_assemblies.json"
            input_manifest.parent.mkdir()
            input_manifest.write_text(json.dumps({
                "empty_sample": {
                    "path": "/inputs/empty.fna", "assembly_type": "metagenome",
                    "organism": "", "sha256": "abc", "size_bytes": 10,
                    "contig_count": 1, "total_length": 1000, "source_layout": "manifest",
                }
            }), encoding="utf-8")
            assembly_dir = root / "amr" / "assemblies" / "empty_sample"
            for _, (input_name, columns) in aggregate_amr.TABLES.items():
                aggregate_amr.write_rows(assembly_dir / input_name, columns, [])
            aggregate_amr.write_rows(
                assembly_dir / "digest.tsv", aggregate_amr.DIGEST_COLUMNS, []
            )
            (assembly_dir / "qc.json").write_text(json.dumps({
                "sources": [
                    {"source": "amrfinderplus", "retained_hits": 0, "hits_without_coordinates": 0},
                    {"source": "rgi", "retained_hits": 0, "hits_without_coordinates": 0},
                    {"source": "genomad", "retained_regions": 0},
                ],
                "loci": 0, "multi_tool_loci": 0, "mobility_links": 0, "mobile_loci": 0,
            }), encoding="utf-8")
            args = argparse.Namespace(
                input_manifest=str(input_manifest),
                per_assembly_dir=str(root / "amr" / "assemblies"),
                output_dir=str(root / "amr"), minimum_overlap=0.8,
                rgi_alignment_tool="DIAMOND", rgi_include_loose=False,
                rgi_include_nudge=False, genomad_preset="default", genomad_splits=8,
                amrfinder_db="/db/amrfinder", card_db="/db/card", genomad_db="/db/genomad",
                drakkar_version="2.3.0",
                amrfinder_version="4.2.7", rgi_version="6.0.8", genomad_version="1.11.0",
                prodigal_version="2.6.3",
            )
            aggregate_amr.aggregate(args)

            with (root / "amr" / "assembly_summary.tsv").open(encoding="utf-8") as handle:
                summaries = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(summaries), 1)
            self.assertEqual(summaries[0]["assembly_id"], "empty_sample")
            self.assertEqual(summaries[0]["amr_loci"], "0")
            with lzma.open(root / "amr" / "amr_hits.tsv.xz", "rt", encoding="utf-8") as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 0)
            manifest = yaml.safe_load((root / "amr" / "manifest.yaml").read_text(encoding="utf-8"))
            self.assertEqual(manifest["schema_version"], "drakkar-amr-manifest-v1")
            self.assertEqual(manifest["outputs"]["amr_hits.tsv.xz"]["rows"], 0)


if __name__ == "__main__":
    unittest.main()
