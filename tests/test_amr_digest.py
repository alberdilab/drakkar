import argparse
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "drakkar" / "workflow" / "scripts" / "amr_digest.py"
SPEC = importlib.util.spec_from_file_location("amr_digest", SCRIPT)
amr_digest = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(amr_digest)


def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class AmrDigestTests(unittest.TestCase):
    def test_overlapping_callers_are_reconciled_and_get_plasmid_context(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            amrfinder = root / "amrfinder.tsv"
            rgi = root / "rgi.txt"
            plasmids = root / "plasmid_summary.tsv"
            plasmid_genes = root / "plasmid_genes.tsv"
            viruses = root / "virus_summary.tsv"
            virus_genes = root / "virus_genes.tsv"
            write_tsv(amrfinder, [{
                "Protein identifier": "contig_1_1", "Contig id": "contig_1",
                "Start": "100", "Stop": "900", "Strand": "+",
                "Gene symbol": "blaTEM-1", "Sequence name": "TEM beta-lactamase",
                "Element type": "AMR", "Element subtype": "AMR", "Class": "BETA-LACTAM",
                "Subclass": "CEPHALOSPORIN", "Method": "ALLELE",
                "% Identity to reference sequence": "100", "% Coverage of reference sequence": "100",
                "Hierarchy node": "NODE:1",
            }], [
                "Protein identifier", "Contig id", "Start", "Stop", "Strand",
                "Gene symbol", "Sequence name", "Element type", "Element subtype", "Class",
                "Subclass", "Method", "% Identity to reference sequence",
                "% Coverage of reference sequence", "Hierarchy node",
            ])
            write_tsv(rgi, [{
                "ORF_ID": "contig_1_1", "Contig": "contig_1", "Start": "110",
                "Stop": "890", "Orientation": "+", "Best_Hit_ARO": "blaTEM-1",
                "ARO": "3000014", "Drug Class": "beta-lactam",
                "Antibiotic": "cephalosporin", "Resistance Mechanism": "antibiotic inactivation",
                "AMR Gene Family": "TEM beta-lactamase", "Cut_Off": "Perfect",
            }], [
                "ORF_ID", "Contig", "Start", "Stop", "Orientation", "Best_Hit_ARO",
                "ARO", "Drug Class", "Antibiotic", "Resistance Mechanism",
                "AMR Gene Family", "Cut_Off",
            ])
            write_tsv(plasmids, [{
                "seq_name": "contig_1", "length": "2000", "topology": "Circular",
                "plasmid_score": "0.99", "fdr": "0.01", "n_genes": "3",
            }], ["seq_name", "length", "topology", "plasmid_score", "fdr", "n_genes"])
            write_tsv(plasmid_genes, [{
                "gene": "contig_1_1", "annotation_amr": "ARO:3000014",
            }], ["gene", "annotation_amr"])
            write_tsv(viruses, [], ["seq_name", "length", "topology", "coordinates", "virus_score"])
            write_tsv(virus_genes, [], ["gene", "annotation_amr"])

            args = argparse.Namespace(
                assembly_id="sample", amrfinder=str(amrfinder), rgi=str(rgi),
                plasmid_summary=str(plasmids), plasmid_genes=str(plasmid_genes),
                virus_summary=str(viruses), virus_genes=str(virus_genes),
                minimum_overlap=0.8,
                hits_output=str(root / "out" / "amr_hits.tsv"),
                loci_output=str(root / "out" / "amr_loci.tsv"),
                regions_output=str(root / "out" / "mobility_regions.tsv"),
                mobility_output=str(root / "out" / "amr_mobility.tsv"),
                drugs_output=str(root / "out" / "amr_drug_classes.tsv"),
                digest_output=str(root / "out" / "digest.tsv"),
                qc_output=str(root / "out" / "qc.json"),
            )
            amr_digest.digest(args)

            hits = read_tsv(root / "out" / "amr_hits.tsv")
            loci = read_tsv(root / "out" / "amr_loci.tsv")
            digest = read_tsv(root / "out" / "digest.tsv")
            regions = read_tsv(root / "out" / "mobility_regions.tsv")
            self.assertEqual(len(hits), 2)
            self.assertEqual({row["source"] for row in hits}, {"amrfinderplus", "rgi"})
            self.assertEqual({row["locus_id"] for row in hits}, {loci[0]["locus_id"]})
            self.assertEqual(loci[0]["support_status"], "amrfinder_and_rgi")
            self.assertEqual(loci[0]["concordance"], "gene_match")
            self.assertEqual(digest[0]["mobility_status"], "plasmid")
            self.assertEqual(regions[0]["amr_accessions"], "ARO:3000014")
            self.assertIn("Protein identifier", json.loads(hits[0]["raw_details"]))

    def test_drug_explosion_does_not_invent_class_subclass_pairs(self):
        hit = {
            "assembly_id": "sample", "hit_id": "hit:1", "source": "rgi",
            "drug_class": "class A;class B", "drug_subclass": "drug 1;drug 2",
            "resistance_mechanism": "efflux", "gene_family": "family", "ontology_id": "ARO:1",
        }
        rows = amr_digest.drug_rows([hit], {"hit:1": "locus:1"})
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["drug_subclass"] for row in rows}, {"drug 1;drug 2"})

    def test_same_source_overlaps_remain_distinct_and_provirus_coordinates_are_mapped(self):
        hits = [
            {
                "assembly_id": "sample", "hit_id": f"hit:{index}",
                "source": "amrfinderplus", "contig": "contig_2", "start": start,
                "end": end, "strand": "+", "gene_symbol": symbol,
                "gene_family": symbol, "ontology_id": "", "drug_class": "",
                "drug_subclass": "", "resistance_mechanism": "",
            }
            for index, (start, end, symbol) in enumerate(
                [(100, 300, "geneA"), (110, 290, "geneB")], start=1
            )
        ]
        loci, _ = amr_digest.reconcile_hits(hits, "sample", 0.8)
        self.assertEqual(len(loci), 2)

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            plasmids = root / "plasmids.tsv"
            plasmid_genes = root / "plasmid_genes.tsv"
            viruses = root / "viruses.tsv"
            virus_genes = root / "virus_genes.tsv"
            write_tsv(plasmids, [], ["seq_name", "length", "plasmid_score"])
            write_tsv(plasmid_genes, [], ["gene", "annotation_amr"])
            write_tsv(viruses, [{
                "seq_name": "contig_2|provirus_1", "length": "401",
                "topology": "Provirus", "coordinates": "500-900", "virus_score": "0.95",
            }], ["seq_name", "length", "topology", "coordinates", "virus_score"])
            write_tsv(virus_genes, [], ["gene", "annotation_amr"])
            regions = amr_digest.parse_genomad(
                plasmids, plasmid_genes, viruses, virus_genes, "sample"
            )
            self.assertEqual(regions[0]["context_type"], "provirus")
            self.assertEqual(regions[0]["contig"], "contig_2")
            self.assertEqual((regions[0]["start"], regions[0]["end"]), (500, 900))


if __name__ == "__main__":
    unittest.main()
