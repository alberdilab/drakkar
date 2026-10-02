from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import yaml


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


def ncbifam_domtbl_row(
    model: str,
    accession: str,
    gene: str,
    *,
    model_length: int,
    query_length: int,
    full_evalue: str,
    full_bitscore: float,
    domain_bitscore: float,
    model_start: int = 2,
    model_end: int = 101,
    query_start: int = 5,
    query_end: int = 104,
    description: str = "model description",
) -> str:
    return " ".join(map(str, [
        model, accession, model_length, gene, "-", query_length,
        full_evalue, full_bitscore, 0.0, 1, 1, full_evalue, full_evalue,
        domain_bitscore, 0.0, model_start, model_end, query_start, query_end,
        query_start, query_end, 0.98, description,
    ]))


def write_ncbifam_metadata(path: Path, *, missing_tigr_cutoff: bool = False) -> Path:
    header = (
        "#ncbi_accession\tsource_identifier\tlabel\tsequence_cutoff\t"
        "domain_cutoff\thmm_length\tfamily_type\tfor_structural_annotation\t"
        "for_naming\tfor_AMRFinder\tproduct_name\tgene_symbol\tgene_synonyms\t"
        "ec_numbers\tgo_terms\tpmids\ttaxonomic_range\t"
        "taxonomic_range_name\ttaxonomic_rank_name\tn_refseq_protein_hits\t"
        "source\tname_orig\thmm_name\tcomment\n"
    )
    tigr_cutoff = "" if missing_tigr_cutoff else "500"
    path.write_text(
        header
        + "NF040708.3\t\tSiroheme_Dcarb_AhbA\t140\t140\t144\tequivalog\tY\tY\tN\t"
        "siroheme decarboxylase subunit alpha\tahbA\t\t4.1.1.111\tGO:0006783\t"
        "21197080\t131567\tcellular organisms\tcellular root\t1773\tNCBIFAM\t\t"
        "siroheme decarboxylase subunit alpha\talpha family\n"
        + f"TIGR04545.1\tTIGR04545\trSAM_ahbD_hemeb\t{tigr_cutoff}\t{tigr_cutoff}\t"
        "339\tequivalog\tY\tY\tN\theme b synthase\tahbD\t\t1.3.98.6\t"
        "GO:0006785\t24713144\t131567\tcellular organisms\tcellular root\t622\t"
        "JCVI\theme b synthase\tAdoMet-dependent heme b synthase\theeme family\n",
        encoding="utf-8",
    )
    return path


class MergeGeneAnnotationTests(unittest.TestCase):
    def test_default_identity_threshold_is_50(self) -> None:
        module = load_merge_module()
        self.assertEqual(module.DEFAULT_IDENTITY_THRESHOLD, 50.0)
        self.assertEqual(module.DEFAULT_QUERY_COVERAGE_THRESHOLD, 0.5)
        self.assertEqual(module.DEFAULT_TARGET_COVERAGE_THRESHOLD, 0.5)
        self.assertNotIn("ncbifam", module.normalize_enabled_sources(None))

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

    def test_ncbifam_preserves_exact_nf_and_versioned_tigr_accessions(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            release = Path(tmpdir) / "ncbifam" / "20.0"
            release.mkdir(parents=True)
            metadata = write_ncbifam_metadata(release / "hmm_PGAP.tsv")
            (release / "database_versions.yaml").write_text(
                yaml.safe_dump({
                    "requested_version": "20.0",
                    "source_version": "NCBIfam/PGAP HMM release 20.0",
                    "sources": ["https://ftp.ncbi.nlm.nih.gov/hmm/20.0/hmm_PGAP.LIB"],
                    "files": [
                        {"filename": "hmm_PGAP.LIB", "sha256": "library-sha"},
                        {"filename": "hmm_PGAP.tsv", "sha256": "metadata-sha"},
                    ],
                }),
                encoding="utf-8",
            )
            raw = Path(tmpdir) / "MAG_A.tblout"
            raw.write_text(
                "# hmmscan :: search sequence(s) against a profile database\n"
                + ncbifam_domtbl_row(
                    "Siroheme_Dcarb_AhbA", "NF040708.3", "c1_1",
                    model_length=144, query_length=200, full_evalue="1e-40",
                    full_bitscore=220, domain_bitscore=210,
                )
                + "\n"
                + ncbifam_domtbl_row(
                    "rSAM_ahbD_hemeb", "TIGR04545.1", "c1_1",
                    model_length=339, query_length=400, full_evalue="1e-80",
                    full_bitscore=560, domain_bitscore=550,
                )
                + "\n",
                encoding="utf-8",
            )

            parsed = module.parse_ncbifam(raw, metadata)

        self.assertEqual(parsed["annotation_id"].tolist(), ["NF040708.3", "TIGR04545.1"])
        self.assertEqual(parsed["source"].tolist(), ["ncbifam", "ncbifam"])
        self.assertEqual(parsed["method"].tolist(), ["hmmer", "hmmer"])
        self.assertEqual(parsed["gene"].tolist(), ["c1_1", "c1_1"])
        self.assertEqual(parsed["hit_rank"].tolist(), [1, 2])
        self.assertEqual(parsed["is_primary"].tolist(), [True, False])
        self.assertNotIn("TIGRFAM", parsed["annotation_id"].tolist())
        details = json.loads(parsed.iloc[1]["details"])
        self.assertEqual(details["source_release"], "20.0")
        self.assertEqual(details["family_type"], "equivalog")
        self.assertEqual(details["profile_source"], "JCVI")
        self.assertEqual(details["native_hmm"]["accession"], "TIGR04545.1")
        self.assertEqual(details["threshold_type"], "trusted_cutoff")
        self.assertEqual(parsed.attrs["annotation_qc"]["database_release"], "20.0")
        self.assertEqual(
            json.loads(parsed.attrs["annotation_qc"]["database_checksums"]),
            {"hmm_PGAP.LIB": "library-sha", "hmm_PGAP.tsv": "metadata-sha"},
        )

        projected = [
            (row.gene, "NCBIFAM", row.annotation_id)
            for row in parsed.itertuples(index=False)
        ]
        self.assertEqual(projected, [
            ("c1_1", "NCBIFAM", "NF040708.3"),
            ("c1_1", "NCBIFAM", "TIGR04545.1"),
        ])

    def test_ncbifam_rechecks_trusted_cutoffs_without_evalue_fallback(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            metadata = write_ncbifam_metadata(tmp / "hmm_PGAP.tsv")
            raw = tmp / "hits.tblout"
            raw.write_text(
                ncbifam_domtbl_row(
                    "Siroheme_Dcarb_AhbA", "NF040708.3", "accepted",
                    model_length=144, query_length=160, full_evalue="1e-20",
                    full_bitscore=150, domain_bitscore=145,
                )
                + "\n"
                + ncbifam_domtbl_row(
                    "rSAM_ahbD_hemeb", "TIGR04545.1", "below_tc",
                    model_length=339, query_length=350, full_evalue="1e-100",
                    full_bitscore=499, domain_bitscore=499,
                )
                + "\n",
                encoding="utf-8",
            )

            parsed = module.parse_ncbifam(raw, metadata)

        self.assertEqual(parsed["gene"].tolist(), ["accepted"])
        self.assertEqual(parsed["threshold"].tolist(), [140])
        self.assertEqual(parsed.attrs["annotation_qc"]["reported_records"], 2)
        self.assertEqual(parsed.attrs["annotation_qc"]["rejected_records"], 1)

    def test_ncbifam_profile_without_trusted_cutoff_is_a_hard_error(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            metadata = write_ncbifam_metadata(
                tmp / "hmm_PGAP.tsv", missing_tigr_cutoff=True
            )
            raw = tmp / "hits.tblout"
            raw.write_text(
                ncbifam_domtbl_row(
                    "rSAM_ahbD_hemeb", "TIGR04545.1", "gene1",
                    model_length=339, query_length=350, full_evalue="1e-100",
                    full_bitscore=550, domain_bitscore=540,
                )
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError, "no complete trusted cutoff.*does not use an E-value fallback"
            ):
                module.parse_ncbifam(raw, metadata)

    def test_ncbifam_does_not_coerce_an_unversioned_tigr_accession(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            metadata = write_ncbifam_metadata(tmp / "hmm_PGAP.tsv")
            raw = tmp / "hits.tblout"
            raw.write_text(
                ncbifam_domtbl_row(
                    "rSAM_ahbD_hemeb", "TIGR04545", "gene1",
                    model_length=339, query_length=350, full_evalue="1e-100",
                    full_bitscore=550, domain_bitscore=540,
                )
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError, "absent from hmm_PGAP.tsv.*TIGR04545"
            ):
                module.parse_ncbifam(raw, metadata)

    def test_ncbifam_rows_survive_the_full_gene_table_merge_losslessly(self) -> None:
        module = load_merge_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            gff = tmp / "genes.gff"
            gff.write_text(
                "c1\tProdigal\tCDS\t1\t1200\t.\t+\t0\tID=1_1;partial=00\n",
                encoding="utf-8",
            )
            metadata = write_ncbifam_metadata(tmp / "hmm_PGAP.tsv")
            raw = tmp / "hits.tblout"
            raw.write_text(
                ncbifam_domtbl_row(
                    "Siroheme_Dcarb_AhbA", "NF040708.3", "c1_1",
                    model_length=144, query_length=400, full_evalue="1e-40",
                    full_bitscore=220, domain_bitscore=210,
                )
                + "\n"
                + ncbifam_domtbl_row(
                    "rSAM_ahbD_hemeb", "TIGR04545.1", "c1_1",
                    model_length=339, query_length=400, full_evalue="1e-80",
                    full_bitscore=560, domain_bitscore=550,
                )
                + "\n",
                encoding="utf-8",
            )
            output = tmp / "genes.tsv"
            qc = tmp / "genes.qc.json"

            merged = module.merge_annotations(
                str(gff), "", "", "", "", "", "", "", "", "", "",
                str(output),
                mag="MAG_A",
                enabled_sources={"ncbifam"},
                qc_output=qc,
                ncbifam_file=raw,
                ncbifam_metadata_file=metadata,
            )
            qc_payload = json.loads(qc.read_text(encoding="utf-8"))
            expected_release = tmp.name

        rows = merged[merged["source"] == "ncbifam"]
        self.assertEqual(rows["annotation_id"].tolist(), ["NF040708.3", "TIGR04545.1"])
        self.assertEqual(
            [(row.gene, "NCBIFAM", row.annotation_id) for row in rows.itertuples()],
            [
                ("c1_1", "NCBIFAM", "NF040708.3"),
                ("c1_1", "NCBIFAM", "TIGR04545.1"),
            ],
        )
        self.assertEqual(
            next(
                record for record in qc_payload["sources"]
                if record["source"] == "ncbifam"
            )["database_release"],
            expected_release,
        )

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

    def test_amr_parser_reads_amrfinderplus_and_ranks_by_its_own_evidence_tier(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            amr = Path(tmpdir) / "amr.tsv"
            amr.write_text(
                "Protein identifier\tContig id\tStart\tStop\tStrand\tElement symbol\t"
                "Element name\tType\tSubtype\tMethod\tClass\tSubclass\t"
                "% Identity to reference\t% Coverage of reference\tHierarchy node\n"
                # A weaker HMM call listed first, to prove ranking is not file order.
                "gene1\tc1\t1\t900\t+\tblaOXA\tclass D beta-lactamase\tAMR\tAMR\t"
                "HMM\tBETA-LACTAM\tCEPHALOSPORIN\t\t\tblaOXA\n"
                "gene1\tc1\t1\t900\t+\tblaOXA-48\tcarbapenemase OXA-48\tAMR\tAMR\t"
                "EXACTP\tBETA-LACTAM\tCARBAPENEM\t100.00\t100.00\tblaOXA-48\n"
                # A "plus" stress gene, which this source deliberately drops.
                "gene2\tc1\t1000\t1600\t-\tarsR\tarsenic repressor\tSTRESS\tMETAL\t"
                "BLASTP\tARSENIC\t\t95.00\t99.00\tarsR\n",
                encoding="utf-8",
            )
            parsed = module.parse_amr(amr)

        self.assertEqual(parsed["annotation_id"].tolist(), ["blaOXA-48", "blaOXA"])
        self.assertEqual(parsed["hit_rank"].tolist(), [1, 2])
        self.assertEqual(parsed["is_primary"].tolist(), [True, False])
        self.assertEqual(parsed["method"].tolist(), ["amrfinderplus", "amrfinderplus"])
        self.assertEqual(parsed["annotation_type"].tolist(), ["BETA-LACTAM", "BETA-LACTAM"])
        primary = parsed.iloc[0]
        self.assertEqual(primary["identity"], 100.0)
        # Coverage is stored as a fraction even though AMRFinderPlus reports a percentage.
        self.assertEqual(primary["coverage"], 1.0)
        self.assertEqual(json.loads(primary["details"])["method"], "EXACTP")
        self.assertEqual(
            json.loads(primary["details"])["threshold_rule"], "amrfinderplus_curated"
        )
        # The STRESS row is excluded from this source but stays visible in QC.
        self.assertNotIn("arsR", parsed["annotation_id"].tolist())
        qc = parsed.attrs["annotation_qc"]
        self.assertEqual(qc["source"], "ncbi_amrfinder")
        self.assertEqual(qc["reported_records"], 3)
        self.assertEqual(qc["retained_records"], 2)
        self.assertEqual(qc["rejected_records"], 1)
        self.assertEqual(qc["filter_stage"], "upstream_native")

    def test_amr_parser_rejects_output_missing_required_columns(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            amr = Path(tmpdir) / "amr.tsv"
            amr.write_text("Contig id\tStart\tStop\nc1\t1\t900\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing required column"):
                module.parse_amr(amr)

    def test_card_parser_keeps_rgi_cutoff_tier_and_curated_threshold(self) -> None:
        module = load_merge_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            card = Path(tmpdir) / "card.txt"
            card.write_text(
                "ORF_ID\tContig\tCut_Off\tBest_Hit_ARO\tARO\tBest_Identities\t"
                "Best_Hit_bit-score\tPass_bit-score\tDrug Class\tAMR Gene Family\t"
                "Resistance Mechanism\tPercentage Length of Reference Sequence\tModel_type\n"
                "gene1 # 1 # 900 # 1 # ID=1_1\tc1\tStrict\tOXA-48\t3001414\t88.5\t"
                "410\t400\tcephalosporin\tOXA beta-lactamase\tantibiotic inactivation\t"
                "97.5\tprotein homolog model\n"
                "gene1 # 1 # 900 # 1 # ID=1_1\tc1\tPerfect\tOXA-181\t3002000\t100.0\t"
                "520\t400\tcarbapenem\tOXA beta-lactamase\tantibiotic inactivation\t"
                "100.0\tprotein homolog model\n",
                encoding="utf-8",
            )
            parsed = module.parse_card(card)

        # Perfect outranks Strict regardless of file order, and the Prodigal
        # coordinates RGI echoes back are stripped from the gene id.
        self.assertEqual(parsed["gene"].tolist(), ["gene1", "gene1"])
        self.assertEqual(parsed["annotation"].tolist(), ["OXA-181", "OXA-48"])
        self.assertEqual(parsed["annotation_id"].tolist(), ["3002000", "3001414"])
        self.assertEqual(parsed["hit_rank"].tolist(), [1, 2])
        primary = parsed.iloc[0]
        self.assertEqual(primary["source"], "card")
        self.assertEqual(primary["method"], "rgi_main")
        self.assertEqual(primary["bitscore"], 520.0)
        self.assertEqual(primary["threshold"], 400.0)
        self.assertEqual(primary["coverage"], 1.0)
        self.assertEqual(json.loads(primary["details"])["cut_off"], "Perfect")
        self.assertEqual(
            json.loads(primary["details"])["resistance_mechanism"],
            "antibiotic inactivation",
        )
        self.assertEqual(parsed.attrs["annotation_qc"]["filter_stage"], "upstream_native")

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
                str(gff), str(kegg), "", str(cutoffs), "", "", "", "", "", "", "",
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
                    str(gff), "", "", "", "", "", "", "", "", "", "", str(tmp / "out.tsv")
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
                    str(gff), str(kegg), "", str(cutoffs), "", "", "", "", "", "", "",
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
                str(gff), str(stale_kegg), "", "", "", "", "", "", "", "",
                str(signalp), str(output), mag="MAG_A", enabled_sources={"signalp"},
            )

        self.assertEqual(set(result["source"]), {"prodigal", "signalp"})


if __name__ == "__main__":
    unittest.main()
