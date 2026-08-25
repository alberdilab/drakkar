from __future__ import annotations

import importlib.util
import json
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


def hmmer_row(model: str, gene: str, evalue: str, bitscore: str, accession: str = "-") -> str:
    return (
        f"{model:<20} {accession:<10} {gene:<12} -           {evalue}  {bitscore}   0.0   "
        f"{evalue}  {bitscore}   0.0   1.0   1   0   0   1   1   1   1 description"
    )


class MergeGeneAnnotationTests(unittest.TestCase):
    def test_default_identity_threshold_is_50(self) -> None:
        module = load_merge_module()
        self.assertEqual(module.DEFAULT_IDENTITY_THRESHOLD, 50.0)
        self.assertEqual(module.DEFAULT_QUERY_COVERAGE_THRESHOLD, 0.5)
        self.assertEqual(module.DEFAULT_TARGET_COVERAGE_THRESHOLD, 0.5)

    def test_kofam_hits_require_the_native_cutoff_table(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            kegg = Path(tmpdir) / "kegg.tsv"
            kegg.write_text(
                "# columns\n" + hmmer_row("K00001", "gene1", "1e-30", "120") + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "native ko_list cutoff table is missing"):
                module.parse_kegg(kegg, "", "", 1e-10)

    def test_kofam_native_score_cutoff_is_authoritative_over_fallback_evalue(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            kegg = tmp / "kegg.tsv"
            cutoffs = tmp / "ko_list.tsv"
            kegg.write_text(
                "# columns\n"
                + hmmer_row("K00001", "gene_native_pass", "1e-5", "120")
                + "\n"
                + hmmer_row("K00002", "gene_native_fail", "1e-30", "50")
                + "\n",
                encoding="utf-8",
            )
            cutoffs.write_text(
                "knum\tthreshold\tscore_type\n"
                "K00001\t100\tfull\n"
                "K00002\t100\tfull\n",
                encoding="utf-8",
            )

            parsed = module.parse_kegg(kegg, "", cutoffs, 1e-10)

        self.assertEqual(parsed["gene"].tolist(), ["gene_native_pass"])
        self.assertEqual(parsed["annotation_id"].tolist(), ["K00001"])
        self.assertEqual(parsed["score"].tolist(), [120])
        self.assertEqual(parsed["score_type"].tolist(), ["full_bitscore"])
        self.assertEqual(parsed["threshold"].tolist(), [100])

    def test_kofam_primary_hit_uses_margin_above_native_cutoff(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            kegg = tmp / "kegg.tsv"
            cutoffs = tmp / "ko_list.tsv"
            kegg.write_text(
                "# columns\n"
                + hmmer_row("K00001", "gene1", "1e-100", "101")
                + "\n"
                + hmmer_row("K00002", "gene1", "1e-20", "150")
                + "\n",
                encoding="utf-8",
            )
            cutoffs.write_text(
                "knum\tthreshold\tscore_type\n"
                "K00001\t100\tfull\n"
                "K00002\t100\tfull\n",
                encoding="utf-8",
            )

            parsed = module.parse_kegg(kegg, "", cutoffs, 1e-10)

        self.assertEqual(parsed["annotation_id"].tolist(), ["K00002", "K00001"])
        self.assertEqual(parsed["rank_score"].tolist(), [50, 1])
        self.assertEqual(parsed["rank_score_type"].tolist(), [
            "bitscore_above_kofam_cutoff", "bitscore_above_kofam_cutoff"
        ])

    def test_repeated_kegg_hierarchy_nodes_do_not_duplicate_a_kofam_hit(self) -> None:
        module = load_merge_module()

        def hierarchy_branch(name: str) -> dict:
            return {
                "children": [
                    {"children": [{"children": [{"name": name}]}]},
                ]
            }

        hierarchy = {
            "children": [
                hierarchy_branch("K00001 first placement [EC:1.1.1.1]"),
                hierarchy_branch("K00001 second placement [EC:2.2.2.2]"),
                hierarchy_branch("K00001 repeated placement [EC:1.1.1.1]"),
            ]
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            kegg = tmp / "kegg.tsv"
            hierarchy_path = tmp / "ko00001.json"
            cutoffs = tmp / "ko_list.tsv"
            kegg.write_text(
                "# columns\n" + hmmer_row("K00001", "gene1", "1e-30", "120") + "\n",
                encoding="utf-8",
            )
            hierarchy_path.write_text(json.dumps(hierarchy), encoding="utf-8")
            cutoffs.write_text(
                "knum\tthreshold\tscore_type\nK00001\t100\tfull\n",
                encoding="utf-8",
            )

            parsed = module.parse_kegg(kegg, hierarchy_path, cutoffs, 1e-10)

        self.assertEqual(len(parsed), 1)
        self.assertEqual(parsed["annotation_id"].tolist(), ["K00001"])
        self.assertEqual(parsed["annotation"].tolist(), ["1.1.1.1;2.2.2.2"])
        self.assertEqual(parsed["hit_rank"].tolist(), [1])
        self.assertEqual(json.loads(parsed.iloc[0]["details"])["ec"], "1.1.1.1;2.2.2.2")

    def test_vfdb_parser_preserves_every_qualifying_hit_and_native_scores(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vf_hits = tmp / "vfdb.txt"
            vf_map = tmp / "vfdb.tsv"
            vf_hits.write_text(
                "gene1\tentry_low_identity\t95\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200\t120\t130\t0.83\t0.77\n"
                "gene1\tentry_second\t99\t90\t1\t0\t2\t91\t5\t94\t1e-20\t180\t120\t140\t0.75\t0.64\n"
                "gene1\tentry_best\t98\t110\t0\t1\t1\t110\t1\t109\t1e-30\t210\t120\t120\t0.92\t0.91\n"
                "gene1\tentry_low_coverage\t100\t20\t0\t0\t1\t20\t1\t20\t1e-80\t250\t100\t100\t0.20\t0.20\n"
                "gene2\tentry_bad_evalue\t99\t100\t0\t0\t1\t100\t1\t100\t1e-5\t120\t100\t100\t1.0\t1.0\n",
                encoding="utf-8",
            )
            vf_map.write_text(
                "entry\tvf\tvfc\tvf_type\tmapping_schema\n"
                "entry_low_identity\tlow\tVFC0001\tlow_type\tdrakkar-vfdb-v2\n"
                "entry_second\tsecond\tVFC0002\tsecond_type\tdrakkar-vfdb-v2\n"
                "entry_best\tbest\tVFC0003\tbest_type\tdrakkar-vfdb-v2\n"
                "entry_low_coverage\tshort\tVFC0005\tshort_type\tdrakkar-vfdb-v2\n"
                "entry_bad_evalue\tbad\tVFC0004\tbad_type\tdrakkar-vfdb-v2\n",
                encoding="utf-8",
            )

            parsed = module.parse_vfdb(vf_hits, vf_map, 1e-10, 98)

        self.assertEqual(parsed["annotation_id"].tolist(), ["entry_best", "entry_second"])
        self.assertEqual(parsed["hit_rank"].tolist(), [1, 2])
        self.assertEqual(parsed["is_primary"].tolist(), [True, False])
        self.assertEqual(parsed["source"].tolist(), ["vfdb", "vfdb"])
        self.assertEqual(parsed["method"].tolist(), ["mmseqs_easy_search", "mmseqs_easy_search"])
        self.assertEqual(parsed["identity"].tolist(), [98, 99])
        self.assertEqual(parsed["bitscore"].tolist(), [210, 180])
        self.assertEqual(parsed["coverage"].tolist(), [0.91, 0.64])
        self.assertEqual(parsed["rank_score_type"].tolist(), [
            "minimum_query_target_coverage", "minimum_query_target_coverage"
        ])
        self.assertEqual(json.loads(parsed.iloc[0]["details"])["vfc"], "VFC0003")
        self.assertEqual(parsed.attrs["annotation_qc"]["rejected_records"], 3)

    def test_vfdb_parser_rejects_fractional_identity(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vf_hits = tmp / "vfdb.txt"
            vf_map = tmp / "vfdb.tsv"
            vf_hits.write_text(
                "gene1\tentry_a\t0.99\t100\t1\t0\t1\t100\t1\t100\t1e-50\t200\t100\t100\t1.0\t1.0\n"
                "gene2\tentry_b\t0.75\t100\t25\t0\t1\t100\t1\t100\t1e-20\t150\t100\t100\t1.0\t1.0\n",
                encoding="utf-8",
            )
            vf_map.write_text("entry\tvf\tvf_type\nentry_a\tA\ttype_a\nentry_b\tB\ttype_b\n", encoding="utf-8")

            parsed = module.parse_vfdb(vf_hits, vf_map, 1e-10, 50.0)

        self.assertEqual(len(parsed), 0, "Fractional identity values must not pass the percentage threshold")

    def test_vfdb_parser_rejects_legacy_mapping_schema(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vf_hits = tmp / "vfdb.txt"
            vf_map = tmp / "vfdb.tsv"
            vf_hits.write_text(
                "gene1\tentry_a\t99\t100\t1\t0\t1\t100\t1\t100\t1e-50\t200\t100\t100\t1.0\t1.0\n",
                encoding="utf-8",
            )
            vf_map.write_text(
                "entry\tvf\tvfc\tvf_type\n"
                "entry_a\tadhesin [Adherence (VFC0001)] [Escherichia coli]"
                "\tVFC0001\tEscherichia coli\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "incompatible with Drakkar 2.0"):
                module.parse_vfdb(vf_hits, vf_map, 1e-10, 50.0)

    def test_cazy_parser_preserves_domains_scores_and_coordinates(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            dbcan = Path(tmpdir) / "dbCAN_hmm_results.tsv"
            dbcan.write_text(
                "HMM Name\tHMM Length\tTarget Name\tTarget Length\ti-Evalue\t"
                "HMM From\tHMM To\tTarget From\tTarget To\tCoverage\tHMM File Name\n"
                "GH5.hmm\t300\tgene1\t500\t1e-20\t1\t150\t20\t170\t0.50\tdbCAN.hmm\n"
                "CBM6.hmm\t120\tgene1\t500\t1e-18\t1\t80\t220\t299\t0.67\tdbCAN.hmm\n"
                "GH5.hmm\t300\tgene1\t500\t1e-25\t1\t150\t320\t469\t0.50\tdbCAN.hmm\n"
                "GT2.hmm\t250\tgene2\t400\t1e-30\t1\t200\t50\t249\t0.80\tdbCAN.hmm\n",
                encoding="utf-8",
            )

            parsed = module.parse_cazy(dbcan)

        gene1 = parsed[parsed["gene"] == "gene1"]
        self.assertEqual(len(gene1), 3)
        self.assertEqual(gene1["annotation_id"].tolist(), ["CBM6", "GH5", "GH5"])
        self.assertEqual(gene1["hit_rank"].tolist(), [1, 2, 3])
        self.assertEqual(gene1["coverage"].tolist(), [0.67, 0.50, 0.50])
        self.assertEqual(gene1["query_start"].tolist(), [220, 320, 20])
        self.assertTrue(all(gene1["evidence"] == "sequence_homology"))

    def test_cazy_parser_rejects_legacy_hmmer_table(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            legacy = Path(tmpdir) / "cazy.tsv"
            legacy.write_text("# hmmscan tblout has no dbCAN coverage column\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing required columns"):
                module.parse_cazy(legacy)

    def test_pfam_parser_preserves_multiple_families_and_all_ec_mappings(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            pfam = tmp / "pfam.tsv"
            ec_map = tmp / "pfam_ec.tsv"
            pfam.write_text(
                "# columns\n"
                + hmmer_row("FamilyA", "gene1", "1e-20", "80", "PF00001.2")
                + "\n"
                + hmmer_row("FamilyB", "gene1", "1e-10", "60", "PF00002.1")
                + "\n",
                encoding="utf-8",
            )
            ec_map.write_text(
                "Type\tConfidence-Score\tPfam-Domain\tEC-Number\n"
                "GOLD\t0.9\tPF00001\t1.1.1.1\n"
                "GOLD\t0.8\tPF00001\t2.2.2.2\n",
                encoding="utf-8",
            )

            parsed = module.parse_pfam(pfam, ec_map)

        self.assertEqual(parsed["annotation_id"].tolist(), ["PF00001", "PF00002"])
        self.assertEqual(parsed["hit_rank"].tolist(), [1, 2])
        associations = json.loads(parsed.iloc[0]["details"])["ec_associations"]
        self.assertEqual([entry["ec"] for entry in associations], ["1.1.1.1", "2.2.2.2"])

    def test_amr_parser_preserves_multiple_models_and_hash_prefixed_mapping_header(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            amr = tmp / "amr.tsv"
            amr_map = tmp / "amr_map.tsv"
            amr.write_text(
                "# columns\n"
                + hmmer_row("ModelA", "gene1", "1e-30", "120", "ACC_A")
                + "\n"
                + hmmer_row("ModelB", "gene1", "1e-15", "90", "ACC_B")
                + "\n",
                encoding="utf-8",
            )
            amr_map.write_text(
                "#hmm_accession\tgene_symbol\tsubtype\tsubclass\n"
                "ACC_A\tblaA\tbeta-lactamase\tbeta-lactam\n"
                "ACC_B\ttetB\tefflux\ttetracycline\n",
                encoding="utf-8",
            )

            parsed = module.parse_amr(amr, amr_map)

        self.assertEqual(parsed["annotation_id"].tolist(), ["ACC_A", "ACC_B"])
        self.assertEqual(parsed["annotation"].tolist(), ["blaA", "tetB"])
        self.assertEqual(parsed["hit_rank"].tolist(), [1, 2])
        self.assertEqual(json.loads(parsed.iloc[1]["details"])["mappings"][0]["subclass"], "tetracycline")

    def test_signalp_and_defensefinder_keep_multiple_predictions(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            signalp = tmp / "signalp.tsv"
            defense = tmp / "defense.tsv"
            signalp.write_text("gene1\tSP\t0.7\ngene1\tLIPO\t0.9\n", encoding="utf-8")
            defense.write_text(
                "hit_id\tgene_name\ttype\tactivity\thit_score\n"
                "gene1\tdefA\tSystemA\tDefense\t50\n"
                "gene1\tantiA\tSystemB\tAntidefense\t40\n",
                encoding="utf-8",
            )

            signal_hits = module.parse_signalp(signalp)
            defense_hits = module.parse_defensefinder(defense)

        self.assertEqual(signal_hits["annotation_id"].tolist(), ["LIPO", "SP"])
        self.assertEqual(signal_hits["hit_rank"].tolist(), [1, 2])
        self.assertEqual(defense_hits["annotation_id"].tolist(), ["defA", "antiA"])
        self.assertEqual(defense_hits["annotation_type"].tolist(), ["Defense", "Antidefense"])

    def test_uniprot_accession_from_target(self) -> None:
        module = load_merge_module()
        self.assertEqual(module.uniprot_accession_from_target("AF-P12345-F1-model_v4.cif.gz"), "P12345")
        self.assertEqual(module.uniprot_accession_from_target("AF-A0A0B1-F1-model_v4"), "A0A0B1")
        self.assertEqual(module.uniprot_accession_from_target("1abc_A"), "1abc_A")

    def test_foldseek_parser_preserves_all_passing_hits_and_mappings(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            m8 = tmp / "foldseek.m8"
            mapping = tmp / "map.tsv"
            m8.write_text(
                "gene1\tAF-P_worse-F1-model_v4.cif.gz\t0.4\t100\t10\t0\t1\t100\t1\t100\t1e-12\t90\n"
                "gene1\tAF-P_keep-F1-model_v4.cif.gz\t0.6\t100\t5\t0\t1\t100\t1\t100\t1e-30\t150\n"
                "gene2\tAF-P_bad-F1-model_v4.cif.gz\t0.3\t100\t40\t0\t1\t100\t1\t100\t1e-5\t40\n",
                encoding="utf-8",
            )
            mapping.write_text(
                "accession\tkegg\tec\tpfam\n"
                "P_keep\tK00010\t1.1.1.1\tPF00010\n"
                "P_worse\tK99999\t9.9.9.9\tPF99999\n",
                encoding="utf-8",
            )

            parsed = module.parse_foldseek(m8, mapping, 1e-10)

        self.assertEqual(parsed["annotation_id"].tolist(), ["P_keep", "P_worse"])
        self.assertEqual(parsed["hit_rank"].tolist(), [1, 2])
        self.assertEqual(parsed["bitscore"].tolist(), [150, 90])
        self.assertEqual(parsed.iloc[0]["annotation"], "kegg=K00010;ec=1.1.1.1;pfam=PF00010")
        self.assertEqual(parsed.iloc[1]["identity"], 0.4)

    def test_long_table_preserves_mag_coordinates_provenance_and_unannotated_genes(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            gff = tmp / "genes.gff"
            kegg = tmp / "kegg.tsv"
            m8 = tmp / "foldseek.m8"
            mapping = tmp / "map.tsv"
            cutoffs = tmp / "ko_list.tsv"
            out = tmp / "out.tsv"
            qc = tmp / "out.qc.json"
            gff.write_text(
                "##gff-version 3\n"
                "c1\tProdigal\tCDS\t1\t90\t.\t+\t0\tID=1_1;partial=00\n"
                "c1\tProdigal\tCDS\t100\t200\t.\t+\t0\tID=1_2;partial=00\n"
                "c1\tProdigal\tCDS\t300\t450\t.\t-\t0\tID=1_3;partial=00\n",
                encoding="utf-8",
            )
            kegg.write_text(
                "# columns\n"
                + hmmer_row("K00001", "c1_1", "1e-30", "100.0")
                + "\n"
                + hmmer_row("K00002", "c1_1", "1e-20", "90.0")
                + "\n",
                encoding="utf-8",
            )
            m8.write_text(
                "c1_2\tAF-P_str-F1-model_v4.cif.gz\t0.5\t100\t10\t0\t1\t100\t1\t100\t1e-25\t140\n",
                encoding="utf-8",
            )
            mapping.write_text("accession\tkegg\tec\tpfam\nP_str\tK00003\t3.3.3.3\tPF00003\n", encoding="utf-8")
            cutoffs.write_text(
                "knum\tthreshold\tscore_type\n"
                "K00001\t80\tfull\n"
                "K00002\t80\tfull\n",
                encoding="utf-8",
            )

            result = module.merge_annotations(
                str(gff), str(kegg), "", str(cutoffs), "", "", "", "", "", "", "", "",
                str(out), foldseek_file=str(m8), foldseekdb_file=str(mapping), mag="MAG_A",
                enabled_sources={"kegg", "structure"}, qc_output=qc,
            )

            written = out.read_text(encoding="utf-8").splitlines()
            qc_payload = json.loads(qc.read_text(encoding="utf-8"))

        self.assertEqual(written[0].split("\t"), module.OUTPUT_COLUMNS)
        self.assertEqual(set(result["mag"]), {"MAG_A"})
        self.assertEqual(len(result[result["source"] == "prodigal"]), 3)
        self.assertEqual(len(result[(result["gene"] == "c1_1") & (result["source"] == "kegg")]), 2)
        unannotated = result[result["gene"] == "c1_3"]
        self.assertEqual(unannotated["source"].tolist(), ["prodigal"])
        self.assertEqual(unannotated.iloc[0]["contig"], "c1")
        self.assertEqual(unannotated.iloc[0]["start"], 300)
        self.assertEqual(unannotated.iloc[0]["strand"], "-")
        structure = result[result["source"] == "uniprot_swissprot"].iloc[0]
        self.assertEqual(structure["method"], "foldseek_prostt5")
        self.assertEqual(structure["evidence"], "structure_homology")
        self.assertNotIn("kegg", result.columns)
        self.assertEqual(qc_payload["schema_version"], "drakkar-gene-annotation-qc-v1")
        self.assertEqual(
            {record["source"] for record in qc_payload["sources"]},
            {"prodigal", "kegg", "uniprot_swissprot"},
        )

    def test_mag_identity_is_required(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            gff = tmp / "genes.gff"
            gff.write_text("", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "MAG identity is required"):
                module.merge_annotations(
                    str(gff), "", "", "", "", "", "", "", "", "", "", "", str(tmp / "out.tsv")
                )

    def test_duplicate_prodigal_gene_ids_are_rejected(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            gff = Path(tmpdir) / "genes.gff"
            gff.write_text(
                "c1\tProdigal\tCDS\t1\t90\t.\t+\t0\tID=1_1\n"
                "c1\tProdigal\tCDS\t100\t190\t.\t+\t0\tID=1_1\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "duplicate derived gene IDs: c1_1"):
                module.parse_gene_calls(gff)

    def test_functional_hits_without_a_matching_prodigal_gene_are_rejected(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            gff = tmp / "genes.gff"
            kegg = tmp / "kegg.tsv"
            cutoffs = tmp / "ko_list.tsv"
            output = tmp / "out.tsv"
            gff.write_text(
                "c1\tProdigal\tCDS\t1\t90\t.\t+\t0\tID=1_1\n",
                encoding="utf-8",
            )
            kegg.write_text(
                "# columns\n" + hmmer_row("K00001", "foreign_gene", "1e-30", "100") + "\n",
                encoding="utf-8",
            )
            cutoffs.write_text(
                "knum\tthreshold\tscore_type\nK00001\t80\tfull\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "MAG 'MAG_A'.*kegg:foreign_gene.*this MAG's Prodigal protein FASTA",
            ):
                module.merge_annotations(
                    str(gff), str(kegg), "", str(cutoffs), "", "", "", "", "", "", "", "",
                    str(output), mag="MAG_A", enabled_sources={"kegg"},
                )

            self.assertFalse(output.exists())

    def test_disabled_stale_source_files_are_ignored(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            gff = tmp / "genes.gff"
            stale_kegg = tmp / "stale_kegg.tsv"
            signalp = tmp / "signalp.tsv"
            output = tmp / "out.tsv"
            gff.write_text(
                "c1\tProdigal\tCDS\t1\t90\t.\t+\t0\tID=1_1\n",
                encoding="utf-8",
            )
            stale_kegg.write_text(
                "# columns\n" + hmmer_row("K00001", "c1_1", "1e-30", "100") + "\n",
                encoding="utf-8",
            )
            signalp.write_text("c1_1\tSP\t0.9\n", encoding="utf-8")

            result = module.merge_annotations(
                str(gff), str(stale_kegg), "", "", "", "", "", "", "", "", "",
                str(signalp), str(output), mag="MAG_A", enabled_sources={"signalp"},
            )

        self.assertEqual(set(result["source"]), {"prodigal", "signalp"})


if __name__ == "__main__":
    unittest.main()
