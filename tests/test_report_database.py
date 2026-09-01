from __future__ import annotations

import lzma
import shutil
import sqlite3
import tempfile
import textwrap
import unittest
from pathlib import Path

from drakkar.report import command as report_command
from drakkar.report import ingest as report_ingest
from drakkar.report.schema import SCHEMA_VERSION, connect, create_schema
from drakkar.report.sources import SectionError, parse_sections, probe


GENE_COLUMNS = [
    "mag", "gene", "contig", "start", "end", "strand", "source", "method",
    "evidence", "hit_rank", "is_primary", "rank_score", "rank_score_type",
    "annotation_id", "annotation", "annotation_type", "evalue", "bitscore",
    "score", "score_type", "threshold", "identity", "coverage",
    "query_coverage", "target_coverage", "confidence", "alignment_length",
    "query_start", "query_end", "target_start", "target_end", "model_start",
    "model_end", "details",
]


def gene_row(**values):
    row = {column: "" for column in GENE_COLUMNS}
    row["details"] = '{"native": "blob"}'
    row.update({key: str(value) for key, value in values.items()})
    return "\t".join(row[column] for column in GENE_COLUMNS)


def write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(content).lstrip(), encoding="utf-8")


def write_xz(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with lzma.open(path, "wt", encoding="utf-8") as handle:
        handle.write(textwrap.dedent(content).lstrip())


class ReportFixtureMixin:
    def build_output_dir(self, root: Path) -> Path:
        write(root / "preprocessing.tsv", """
            sample\treads_pre_fastp\tbases_pre_fastp\tmetagenomic_reads\tsinglem_fraction\tnonpareil_diversity
            S1\t1000\t150000\t880\t0.85\t18.2
            S2\t2000\t300000\t1960\tNA\t19.1
            """)
        write(root / "cataloging.tsv", """
            assembly\tsamples\tcoverage_samples\tassembly_contigs\tmapping_rate_percent\tsample_mapping_rates\tmetabat2_bins\tmaxbin2_bins\tsemibin2_bins\tcomebin_bins\tfinal_bins
            A1\tS1,S2\tS1,S2\t1200\t90.0\tS1:88.50;S2:91.20\t5\t4\t6\t3\t7
            A2\tS2\tNA\t800\t0\tNA\t3\t2\t4\t2\t4
            """)
        # Binette's per-assembly final bin report, copied verbatim into
        # cataloging/final. `origin` names the binner a bin came from; a bin
        # several binners produced identically names all of them, and one
        # Binette built itself is marked `binette` and is not original.
        write(root / "cataloging/final/A1.tsv", """
            name\torigin\tis_original\toriginal_name\tcompleteness\tcontamination\tscore\tcheckm2_model\tsize\tN50\tcoding_density\tcontig_count
            binette_bin1\tmetabat\tTrue\tmetabat_3\t95.1\t2.0\t91.1\tgeneral\t3000000\t50000\t0.88\t120
            binette_bin2\tmetabat;semibin\tTrue\tmetabat_5\t88.0\t1.0\t86.0\tgeneral\t2600000\t41000\t0.89\t140
            binette_bin3\tbinette\tFalse\tbinette_bin3\t78.4\t4.5\t69.4\tgeneral\t2100000\t22000\t0.9\t310
            """)
        write(root / "cataloging/final/A2.tsv", """
            name\torigin\tis_original\toriginal_name\tcompleteness\tcontamination\tscore\tcheckm2_model\tsize\tN50\tcoding_density\tcontig_count
            binette_bin1\tmaxbin\tTrue\tmaxbin_1\t91.0\t3.0\t85.0\tgeneral\t2400000\t31000\t0.88\t160
            """)
        write(root / "dereplicating.tsv", """
            input_bin_number\tinput_bin_completeness\tdereplication_ani\toutput_bin_number
            11\t82.4\t0.98\t6
            """)
        # dRep's own data tables. Cdb/Wdb name bins bare, Ndb/Mdb by full path,
        # and Ndb holds both directions of every comparison.
        write(root / "dereplicating/drep/data_tables/Cdb.csv", """
            genome,secondary_cluster,threshold,cluster_method,comparison_algorithm,primary_cluster
            A1_bin_1.fa,1_1,0.02,average,ANImf,1
            A1_bin_2.fa,1_1,0.02,average,ANImf,1
            A1_bin_3.fa,2_1,0.02,average,ANImf,2
            A2_bin_1.fa,2_2,0.02,average,ANImf,2
            A2_bin_2.fa,3_1,0.02,average,ANImf,3
            """)
        write(root / "dereplicating/drep/data_tables/Wdb.csv", """
            genome,cluster,score
            A1_bin_1.fa,1_1,95.0
            A1_bin_3.fa,2_1,91.0
            A2_bin_1.fa,2_2,88.0
            A2_bin_2.fa,3_1,80.0
            """)
        write(root / "dereplicating/drep/data_tables/Ndb.csv", """
            reference,querry,ani,alignment_coverage,primary_cluster
            /scratch/genomes/A1_bin_1.fa,/scratch/genomes/A1_bin_2.fa,0.9990,0.82,1
            /scratch/genomes/A1_bin_2.fa,/scratch/genomes/A1_bin_1.fa,0.9970,0.80,1
            /scratch/genomes/A1_bin_3.fa,/scratch/genomes/A2_bin_1.fa,0.9600,0.71,2
            /scratch/genomes/A2_bin_1.fa,/scratch/genomes/A1_bin_3.fa,0.9600,0.70,2
            """)
        write(root / "dereplicating/drep/data_tables/Mdb.csv", """
            genome1,genome2,dist,similarity
            /scratch/genomes/A1_bin_1.fa,/scratch/genomes/A1_bin_1.fa,0.0,1.0
            /scratch/genomes/A1_bin_1.fa,/scratch/genomes/A1_bin_2.fa,0.02,0.98
            /scratch/genomes/A1_bin_2.fa,/scratch/genomes/A1_bin_1.fa,0.02,0.98
            /scratch/genomes/A1_bin_3.fa,/scratch/genomes/A2_bin_1.fa,0.06,0.94
            /scratch/genomes/A2_bin_1.fa,/scratch/genomes/A1_bin_3.fa,0.06,0.94
            /scratch/genomes/A1_bin_1.fa,/scratch/genomes/A1_bin_3.fa,0.41,0.59
            /scratch/genomes/A1_bin_3.fa,/scratch/genomes/A1_bin_1.fa,0.41,0.59
            /scratch/genomes/A2_bin_2.fa,/scratch/genomes/A1_bin_1.fa,0.45,0.55
            /scratch/genomes/A1_bin_1.fa,/scratch/genomes/A2_bin_2.fa,0.45,0.55
            """)
        write(root / "profiling_genomes/final/mags.tsv", """
            magid\tcompleteness\tcontamination\tsize\tcontigs\tn50\tgc\tcluster\tcluster_members\tscore
            MAG_A\t95.1\t1.5\t2400000\t120\t12000\t44.2\t1_1\t3\t93.6
            MAG_B\tNA\tNA\t2100000\t180\t9000\tNA\t2_1\t2\tNA
            """)
        write(root / "profiling_genomes/final/counts.tsv", """
            Genome\tS1\tS2
            MAG_A\t15000\t22000
            MAG_B\t8000\t5000
            """)
        write(root / "profiling_genomes/final/bases.tsv", """
            Genome\tS1\tS2
            MAG_A\t2100000\t2350000
            MAG_B\t1700000\t1200000
            """)
        # Raw GTDB-Tk concatenation: packed lineage, blob columns, and the
        # archaeal summary's repeated header row.
        write(root / "annotating/genome_taxonomy.tsv", """
            user_genome\tclassification\tclosest_genome_reference\tclosest_genome_ani\tclosest_genome_af\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\tred_value\twarnings
            MAG_A\td__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia;s__Blautia wexlerae\tGCF_000159015.1\t97.8\t0.85\ttopology and ANI\tN/A\tGCF_001.1, Blautia obeum, 95.0, 88.1, 0.6\t92.4\tNA\t
            user_genome\tclassification\tclosest_genome_reference\tclosest_genome_ani\tclosest_genome_af\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\tred_value\twarnings
            MAG_B\td__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__\tGCF_000016525.1\t98.5\t0.91\ttopology and ANI\tN/A\t\t94.0\t0.87\tGenome has warnings
            """)
        write(root / "annotating/annotation_qc.tsv", """
            mag\tlevel\tsource\treported_records\tretained_records\trejected_records\tunmapped_records\tunique_entities\tfilter_stage
            MAG_A\tgene\tkegg\t500\t420\t80\t0\t410\tmerge
            """)
        write_xz(
            root / "annotating/gene_annotations.tsv.xz",
            "\t".join(GENE_COLUMNS) + "\n" + "\n".join([
                gene_row(mag="MAG_A", gene="c1_1", contig="c1", start=1, end=900,
                         strand="+", source="prodigal", hit_rank=1, is_primary="True",
                         annotation_id="CDS"),
                gene_row(mag="MAG_A", gene="c1_1", contig="c1", start=1, end=900,
                         strand="+", source="kegg", hit_rank=1, is_primary="True",
                         annotation_id="K00001", annotation="alcohol dehydrogenase",
                         annotation_type="ko", evalue="1e-40", bitscore="180.2"),
                gene_row(mag="MAG_A", gene="c1_1", contig="c1", start=1, end=900,
                         strand="+", source="kegg", hit_rank=2, is_primary="False",
                         annotation_id="K00002", annotation="aldehyde reductase",
                         annotation_type="ko", evalue="1e-20", bitscore="95.1"),
                # Negative strand must survive: "-" is a value, not a null marker.
                gene_row(mag="MAG_A", gene="c1_2", contig="c1", start=1000, end=1600,
                         strand="-", source="prodigal", hit_rank=1, is_primary="True",
                         annotation_id="CDS"),
                gene_row(mag="MAG_B", gene="c9_4", contig="c9", start=50, end=800,
                         strand="+", source="kegg", hit_rank=1, is_primary="True",
                         annotation_id="K00001", annotation="alcohol dehydrogenase",
                         annotation_type="ko", evalue="1e-35", bitscore="160.0"),
            ]) + "\n",
        )
        write_xz(root / "annotating/cluster_annotations.tsv.xz", """
            mag\tcluster_id\tcontig\tstart\tend\tsource\tmethod\tevidence\ttype\tgene_count\tsubstrate\tgene_functions\tpul_id\tdetails
            MAG_A\tCGC1\tc1\t1\t5000\tdbcan\trun_dbcan_cgc\tcgc_prediction\tCAZyme\t6\tstarch\tGH13,GT2\tPUL0012\t{}
            """)
        write_xz(root / "expressing/gene_counts.tsv.xz", """
            # Program:featureCounts v2.0.6; Command:featureCounts -F GFF
            Geneid\tChr\tStart\tEnd\tStrand\tLength\tS1\tS2
            MAG_A_c1_1\tc1;c1\t1;400\t300;900\t+;+\t800\t120\t340
            MAG_A_c1_2\tc1\t1000\t1600\t-\t600\t45\t0
            """)
        self.write_amr(root)
        write(root / "logging" / "drakkar_20260825-101500.yaml", """
            run_id: '20260825-101500'
            drakkar_version: 2.0.0
            started_at: '2026-08-25T10:15:00+00:00'
            finished_at: '2026-08-25T14:02:11+00:00'
            command: complete
            modules:
            - preprocessing
            - cataloging
            status: completed
            output_directory: /scratch/run
            argv:
            - drakkar
            - complete
            """)
        self.write_benchmark(root, "20260825-101500")
        self.write_failures(root, "20260825-101500")
        return root

    def build_legacy_output_dir(self, root: Path) -> Path:
        """The same directory in the layout used before ``logging/`` existed."""
        self.build_output_dir(root)
        logging_root = root / "logging"
        (logging_root / "drakkar_20260825-101500.yaml").rename(
            root / "drakkar_20260825-101500.yaml"
        )
        shutil.rmtree(logging_root)
        self.write_benchmark(root, "20260825-101500", legacy=True)
        self.write_failures(root, "20260825-101500", legacy=True)
        return root

    def write_amr(self, root: Path) -> None:
        """The tables aggregate_amr writes, in the column order it writes them."""
        write(root / "amr/assembly_summary.tsv", """
            assembly_id\tassembly_type\torganism\tinput_path\tinput_sha256\tinput_size_bytes\tcontig_count\ttotal_length\tamrfinder_hits\trgi_hits\tmobility_regions\tamr_loci\tmulti_tool_loci\tmobile_loci
            A1\tmetagenome\t\t/scratch/A1.fna\tabc123\t5200000\t1200\t5000000\t3\t2\t2\t3\t1\t2
            A2\tisolate\tEscherichia coli\t/scratch/A2.fna\tdef456\t4900000\t80\t4800000\t1\t0\t0\t1\t0\t0
            """)
        write(root / "amr/amr_qc.tsv", """
            assembly_id\tamrfinder_hits\tamrfinder_hits_without_coordinates\trgi_hits\trgi_hits_without_coordinates\tmobility_regions\tamr_loci\tmulti_tool_loci\tmobility_links\tmobile_loci
            A1\t3\t0\t2\t1\t2\t3\t1\t2\t2
            A2\t1\t1\t0\t0\t0\t1\t0\t0\t0
            """)
        write_xz(root / "amr/amr_hits.tsv.xz", """
            assembly_id\thit_id\tlocus_id\tsource\tsource_hit_id\tgene_id\tcontig\tstart\tend\tstrand\tgene_symbol\tgene_name\tontology_id\tdrug_class\tdrug_subclass\tresistance_mechanism\tgene_family\tmodel_type\tdetection_grade\tmethod\tmutation\tidentity\treference_coverage\tbitscore\tthreshold\tis_partial\traw_details
            A1\tH1\tL1\tamrfinderplus\t1\tc1_1\tc1\t100\t900\t+\tblaTEM\tbeta-lactamase TEM\t\tbeta-lactam\t\tantibiotic inactivation\tblaTEM\t\tEXACTX\tEXACTX\t\t100.0\t100.0\t540.0\t\tfalse\t{"native": "blob"}
            A1\tH2\tL1\trgi\t2\tc1_1\tc1\t105\t900\t+\tblaTEM\tTEM-1\tARO:3000873\tbeta-lactam\t\tantibiotic inactivation\tblaTEM\tprotein homolog\tPerfect\tDIAMOND\t\t99.6\t98.0\t530.0\t500\tfalse\t{"native": "blob"}
            A1\tH3\tL2\tamrfinderplus\t3\tc2_4\tc2\t50\t700\t+\tsul1\tsulfonamide resistance\t\tsulfonamide\t\ttarget replacement\tsul1\t\tBLASTX\tBLASTX\t\t96.4\t91.0\t410.0\t\ttrue\t{"native": "blob"}
            A2\tH4\tL4\tamrfinderplus\t1\td1_2\td1\t100\t800\t+\tblaTEM\tbeta-lactamase TEM\t\tbeta-lactam\t\tantibiotic inactivation\tblaTEM\t\tBLASTX\tBLASTX\t\t97.2\t92.0\t480.0\t\tfalse\t{"native": "blob"}
            """)
        write_xz(root / "amr/amr_loci.tsv.xz", """
            assembly_id\tlocus_id\tcontig\tstart\tend\tstrand\tprimary_gene\tgene_symbols\tgene_families\tontology_ids\tdrug_classes\tdrug_subclasses\tresistance_mechanisms\tsources\tsource_count\thit_count\tsupport_status\tconcordance\tdetails
            A1\tL1\tc1\t100\t900\t+\tblaTEM\tblaTEM\tblaTEM\tARO:3000873\tbeta-lactam\t\tantibiotic inactivation\tamrfinderplus;rgi\t2\t2\tamrfinder_and_rgi\tgene_match\t{"hit_ids": ["H1", "H2"]}
            A1\tL2\tc2\t50\t700\t+\tsul1\tsul1\tsul1\t\tsulfonamide\t\ttarget replacement\tamrfinderplus\t1\t1\tamrfinder_only\tsingle_source\t{"hit_ids": ["H3"]}
            A1\tL3\tc3\t10\t400\t-\tqnrS\tqnrS\tqnr\t\tquinolone\t\ttarget protection\trgi\t1\t1\trgi_only\tsingle_source\t{"hit_ids": ["H5"]}
            A2\tL4\td1\t100\t800\t+\tblaTEM\tblaTEM\tblaTEM\t\tbeta-lactam\t\tantibiotic inactivation\tamrfinderplus\t1\t1\tamrfinder_only\tsingle_source\t{"hit_ids": ["H4"]}
            """)
        write_xz(root / "amr/amr_drug_classes.tsv.xz", """
            assembly_id\tlocus_id\thit_id\tsource\tdrug_class\tdrug_subclass\tresistance_mechanism\tgene_family\tontology_id
            A1\tL1\tH1\tamrfinderplus\tbeta-lactam\t\tantibiotic inactivation\tblaTEM\t
            A1\tL1\tH2\trgi\tbeta-lactam\t\tantibiotic inactivation\tblaTEM\tARO:3000873
            A1\tL2\tH3\tamrfinderplus\tsulfonamide\t\ttarget replacement\tsul1\t
            A1\tL3\tH5\trgi\tquinolone\t\ttarget protection\tqnr\t
            A2\tL4\tH4\tamrfinderplus\tbeta-lactam\t\tantibiotic inactivation\tblaTEM\t
            """)
        write_xz(root / "amr/mobility_regions.tsv.xz", """
            assembly_id\tregion_id\tcontig\tstart\tend\tcontext_type\tseq_name\tlength\ttopology\tscore\tfdr\tmarker_enrichment\thallmark_count\tgene_count\tconjugation_genes\tamr_accessions\ttaxonomy
            A1\tR1\tc1\t1\t40000\tplasmid\tc1\t40000\tcircular\t0.94\t0.01\t2.4\t6\t45\ttraB;traC\tARO:3000873\t
            A1\tR2\tc3\t1\t22000\tprovirus\tc3\t22000\tlinear\t0.81\t0.04\t1.2\t3\t28\t\t\tCaudoviricetes
            """)
        write_xz(root / "amr/amr_mobility.tsv.xz", """
            assembly_id\tlocus_id\tregion_id\tcontext_type\tcontig\toverlap_bp\tlocus_overlap_fraction\tregion_score
            A1\tL1\tR1\tplasmid\tc1\t801\t1.0\t0.94
            A1\tL3\tR2\tprovirus\tc3\t391\t1.0\t0.81
            """)

    def write_benchmark(self, root: Path, run_id: str, legacy: bool = False) -> None:
        """The artefacts drakkar.benchmark writes after a SLURM run.

        ``legacy`` writes them where runs predating the logging directory put
        them: the roll-up in the output root, the tables in ``benchmark/``.
        """
        benchmark_dir = root / "benchmark" if legacy else root / "logging" / "benchmark"
        summary_name = (
            f"drakkar_{run_id}_resources.yaml" if legacy else f"drakkar_{run_id}.resources.yaml"
        )
        summary_path = root / summary_name if legacy else benchmark_dir / summary_name
        write(summary_path, f"""
            run_id: '{run_id}'
            command: complete
            profile: slurm
            generated_at: '2026-08-25T14:03:00+00:00'
            status: generated
            benchmarked_launches: 3
            logical_jobs: 2
            retries: 1
            failed_launches: 1
            oom_launches: 1
            timeout_launches: 0
            jobs_missing_accounting: 0
            max_alloc_cpus: 8
            peak_max_rss_mb: 30000.0
            total_elapsed_sec: 10800
            allocated_cpu_sec: 57600
            used_cpu_sec: 43200
            weighted_cpu_efficiency: 0.75
            """)
        write(benchmark_dir / f"drakkar_{run_id}.jobs.tsv", """
            launch_index	rule	attempt	logical_job_key	internal_jobid	external_jobid	wildcards	requested_cpus	requested_mem_mb	requested_runtime_min	state	exit_code	alloc_cpus	elapsed_sec	cpu_time_sec	max_rss_mb	cpu_efficiency	memory_efficiency	runtime_efficiency	oom	timeout
            1	assembly	1	assembly|wildcards|assembly=A1	12	9001	assembly=A1	8	65536	720	OUT_OF_MEMORY	0:125	8	1800	12000	65000.0	0.833	0.9918	0.0417	True	False
            2	assembly	2	assembly|wildcards|assembly=A1	12	9002	assembly=A1	8	131072	720	COMPLETED	0:0	8	5400	38000	98000.0	0.8796	0.7477	0.125	False	False
            3	binning	1	binning|wildcards|assembly=A1	14	9003	assembly=A1	4	16384	240	COMPLETED	0:0	4	3600	9000	8000.0	0.625	0.4883	0.25	False	False
            """)
        write(benchmark_dir / f"drakkar_{run_id}.rules.tsv", """
            rule	launches	logical_jobs	retries	failed_launches	oom_launches	timeout_launches	median_requested_cpus	median_alloc_cpus	median_requested_mem_mb	median_max_rss_mb	median_memory_efficiency	median_requested_runtime_min	median_elapsed_sec	median_runtime_efficiency	allocated_cpu_sec	used_cpu_sec	weighted_cpu_efficiency
            assembly	2	1	1	1	1	0	8	8	98304	81500.0	0.8697	720	3600	0.0833	57600	50000	0.868
            binning	1	1	0	0	0	0	4	4	16384	8000.0	0.4883	240	3600	0.25	14400	9000	0.625
            """)


    def write_failures(self, root: Path, run_id: str, legacy: bool = False) -> None:
        """The failure table drakkar.failures writes after a run that lost a job.

        One row per rule and target: an assembly that ran out of memory twice
        before a retry got it through, a sample whose rule never succeeded, and
        a workflow-level error, which carries ``workflow`` as its target.
        """
        name = (
            f"drakkar_{run_id}_failures.tsv" if legacy
            else f"drakkar_{run_id}.failures.tsv"
        )
        path = root / name if legacy else root / "logging" / name
        write(path, f"""
            run_id	rule	target	attempts	status	category	reason	slurm_state	internal_jobid	external_jobid	detail	action	job_log	output	last_failure_at
            {run_id}	assembly	A1	2	recovered	out-of-memory	the job was killed after exceeding its memory allocation	OUT_OF_MEMORY	12	9001	exceeded its memory allocation (65536 MB requested)	Relaunch the same drakkar command with a larger --memory-multiplier (e.g. --memory-multiplier 2).	logging/slurm/9001.log	cataloging/assembly/A1.fna	2026-08-25T11:40:00+00:00
            {run_id}	singlem	S2	1	failed	timeout	the job hit its SLURM wall-time limit	TIMEOUT	41	9010	hit the SLURM time limit (72 min requested)	Relaunch the same drakkar command with a larger --time-multiplier (e.g. --time-multiplier 2).	logging/slurm/9010.log	preprocessing/singlem/S2_cond.tsv	2026-08-25T12:05:00+00:00
            {run_id}	bowtie2	workflow	1	failed	missing-input	required input files were not found				MissingInputException in rule bowtie2 reference/host.fna.gz	Check the sample information file, input directory, and database paths before relaunching.			
            """)


class SectionParsingTests(unittest.TestCase):
    def test_none_selects_every_section(self):
        self.assertEqual(parse_sections(None)[0], "preprocessing")
        self.assertIn("function", parse_sections(None))

    def test_all_keyword_selects_every_section(self):
        self.assertEqual(parse_sections("all"), parse_sections(None))

    def test_selection_uses_canonical_order(self):
        self.assertEqual(
            parse_sections("function,preprocessing"),
            ("preprocessing", "function"),
        )

    def test_unknown_section_is_rejected_by_name(self):
        with self.assertRaises(SectionError) as raised:
            parse_sections("preprocessing,bogus")
        self.assertIn("bogus", str(raised.exception))

    def test_whitespace_and_case_are_tolerated(self):
        self.assertEqual(parse_sections(" Taxonomy , FUNCTION "), ("taxonomy", "function"))


class ProbeTests(ReportFixtureMixin, unittest.TestCase):
    def test_complete_directory_reports_every_section_available(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            for entry in probe(root):
                self.assertTrue(entry["available"], entry["section"])

    def test_partial_directory_names_the_missing_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root / "preprocessing.tsv", "sample\treads_pre_fastp\nS1\t10\n")
            results = {entry["section"]: entry for entry in probe(root)}
            self.assertTrue(results["preprocessing"]["available"])
            self.assertFalse(results["taxonomy"]["available"])
            self.assertIn(
                "annotating/genome_taxonomy.tsv", results["taxonomy"]["missing"]
            )

    def test_resources_lists_the_benchmark_files_it_found(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            entry = {item["section"]: item for item in probe(root)}["resources"]
            self.assertTrue(entry["available"])
            self.assertIn(
                str(Path("logging") / "drakkar_20260825-101500.yaml"), entry["present"]
            )
            self.assertIn(
                str(Path("logging") / "benchmark" / "drakkar_20260825-101500.resources.yaml"),
                entry["present"],
            )
            self.assertIn(
                str(Path("logging") / "benchmark" / "drakkar_20260825-101500.jobs.tsv"),
                entry["present"],
            )
            self.assertEqual(entry["missing"], [])

    def test_resources_lists_the_failure_table_of_a_run_that_lost_jobs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            entry = {item["section"]: item for item in probe(root)}["resources"]
            self.assertIn(
                str(Path("logging") / "drakkar_20260825-101500.failures.tsv"),
                entry["present"],
            )

    def test_a_run_without_failures_is_not_asked_for_a_failure_table(self):
        # No file is written when nothing failed, so its absence is the good
        # case and must not be reported as something the section is missing.
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            (root / "logging" / "drakkar_20260825-101500.failures.tsv").unlink()
            entry = {item["section"]: item for item in probe(root)}["resources"]
            self.assertTrue(entry["available"])
            self.assertEqual(entry["missing"], [])

    def test_resources_names_the_benchmark_a_run_never_produced(self):
        # A local run has metadata but no benchmark: the section still builds,
        # and the absence is named rather than hidden.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root / "logging" / "drakkar_20260825-101500.yaml", "run_id: '20260825-101500'\n")
            entry = {item["section"]: item for item in probe(root)}["resources"]
            self.assertTrue(entry["available"])
            self.assertEqual(entry["missing"], ["drakkar_<run_id>.resources.yaml"])

    def test_empty_file_does_not_count_as_available(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "preprocessing.tsv").write_text("", encoding="utf-8")
            results = {entry["section"]: entry for entry in probe(root)}
            self.assertFalse(results["preprocessing"]["available"])


class DatabaseBuildTests(ReportFixtureMixin, unittest.TestCase):
    def build(self, root, sections=None, **kwargs):
        db_path = Path(root) / "drakkar.db"
        report_command.build_database(
            root, sections or parse_sections(None), db_path, **kwargs
        )
        connection = sqlite3.connect(db_path)
        connection.row_factory = sqlite3.Row
        return connection

    def test_schema_version_is_stamped_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            rows = connection.execute("SELECT version FROM schema_version").fetchall()
            self.assertEqual([row[0] for row in rows], [SCHEMA_VERSION])
            connection.close()

    def test_gtdbtk_lineage_is_split_into_ranks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT domain, phylum, genus, species, ani, warnings '
                'FROM genome_taxonomy WHERE genome_id = "MAG_A"'
            ).fetchone()
            self.assertEqual(row["domain"], "Bacteria")
            self.assertEqual(row["phylum"], "Bacillota")
            self.assertEqual(row["genus"], "Blautia")
            self.assertEqual(row["species"], "Blautia wexlerae")
            self.assertAlmostEqual(row["ani"], 97.8)
            connection.close()

    def test_repeated_gtdbtk_header_row_is_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            ids = [
                row[0]
                for row in connection.execute("SELECT genome_id FROM genome_taxonomy")
            ]
            self.assertEqual(sorted(ids), ["MAG_A", "MAG_B"])
            connection.close()

    def test_empty_species_rank_becomes_null(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT species FROM genome_taxonomy WHERE genome_id = "MAG_B"'
            ).fetchone()
            self.assertIsNone(row["species"])
            connection.close()

    def test_gene_coordinates_are_stored_once_per_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            genes = connection.execute(
                "SELECT mag, gene, strand FROM gene ORDER BY mag, gene"
            ).fetchall()
            self.assertEqual(
                [(row["mag"], row["gene"]) for row in genes],
                [("MAG_A", "c1_1"), ("MAG_A", "c1_2"), ("MAG_B", "c9_4")],
            )
            connection.close()

    def test_negative_strand_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT strand FROM gene WHERE mag = "MAG_A" AND gene = "c1_2"'
            ).fetchone()
            self.assertEqual(row["strand"], "-")
            connection.close()

    def test_prodigal_rows_are_not_stored_as_annotations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            sources = [
                row[0]
                for row in connection.execute(
                    "SELECT DISTINCT source FROM gene_annotation"
                )
            ]
            self.assertNotIn("prodigal", sources)
            self.assertIn("kegg", sources)
            connection.close()

    def test_annotation_labels_are_deduplicated_into_terms(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            # K00001 is hit by two different genes but stored as one term.
            hits = connection.execute(
                'SELECT COUNT(*) FROM gene_annotation WHERE annotation_id = "K00001"'
            ).fetchone()[0]
            terms = connection.execute(
                'SELECT COUNT(*) FROM annotation_term WHERE annotation_id = "K00001"'
            ).fetchone()[0]
            self.assertEqual(hits, 2)
            self.assertEqual(terms, 1)
            row = connection.execute(
                'SELECT annotation FROM annotation_term WHERE annotation_id = "K00001"'
            ).fetchone()
            self.assertEqual(row["annotation"], "alcohol dehydrogenase")
            connection.close()

    def test_details_column_is_not_carried_into_the_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            columns = {
                row[1]
                for row in connection.execute("PRAGMA table_info(gene_annotation)")
            }
            self.assertNotIn("details", columns)
            connection.close()

    def test_primary_hits_only_drops_secondary_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root, primary_hits_only=True)
            ranks = [
                row[0] for row in connection.execute("SELECT hit_rank FROM gene_annotation")
            ]
            self.assertEqual(set(ranks), {1})
            connection.close()

    def test_secondary_hits_are_kept_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            ranks = sorted(
                row[0] for row in connection.execute("SELECT hit_rank FROM gene_annotation")
            )
            self.assertIn(2, ranks)
            connection.close()

    def test_packed_sample_mapping_rates_are_normalized(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            rows = connection.execute(
                'SELECT sample_id, mapping_rate_percent FROM assembly_sample '
                'WHERE assembly_id = "A1" ORDER BY sample_id'
            ).fetchall()
            self.assertEqual(
                [(row["sample_id"], row["mapping_rate_percent"]) for row in rows],
                [("S1", 88.5), ("S2", 91.2)],
            )
            connection.close()

    def test_binner_counts_are_normalized(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            rows = dict(
                connection.execute(
                    'SELECT binner, bin_count FROM assembly_binner WHERE assembly_id = "A1"'
                ).fetchall()
            )
            self.assertEqual(
                rows, {"metabat2": 5, "maxbin2": 4, "semibin2": 6, "comebin": 3}
            )
            connection.close()

    def test_counts_and_bases_are_melted_into_long_form(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT read_count, covered_bases FROM genome_count '
                'WHERE genome_id = "MAG_A" AND sample_id = "S2"'
            ).fetchone()
            self.assertEqual(row["read_count"], 22000.0)
            self.assertEqual(row["covered_bases"], 2350000.0)
            connection.close()

    def test_na_values_become_null_not_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT completeness, gc FROM genome WHERE genome_id = "MAG_B"'
            ).fetchone()
            self.assertIsNone(row["completeness"])
            self.assertIsNone(row["gc"])
            connection.close()

    def test_featurecounts_comment_line_and_packed_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT contig, start, end, strand, length FROM expressed_gene '
                'WHERE gene_id = "MAG_A_c1_1"'
            ).fetchone()
            # "1;400" / "300;900" collapse to the meta-feature bounds.
            self.assertEqual((row["contig"], row["start"], row["end"]), ("c1", 1, 900))
            self.assertEqual(row["strand"], "+")
            self.assertEqual(row["length"], 800)
            counts = dict(
                connection.execute(
                    'SELECT sample_id, count FROM gene_expression '
                    'WHERE gene_id = "MAG_A_c1_1"'
                ).fetchall()
            )
            self.assertEqual(counts, {"S1": 120.0, "S2": 340.0})
            connection.close()

    def test_run_metadata_is_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute("SELECT * FROM run").fetchone()
            self.assertEqual(row["run_id"], "20260825-101500")
            self.assertEqual(row["command"], "complete")
            self.assertEqual(row["modules"], "preprocessing,cataloging")
            self.assertEqual(row["status"], "completed")
            connection.close()

    def test_benchmark_rollup_is_recorded_per_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute("SELECT * FROM run_benchmark").fetchone()
            self.assertEqual(row["run_id"], "20260825-101500")
            self.assertEqual(row["status"], "generated")
            self.assertEqual(row["profile"], "slurm")
            self.assertEqual(row["benchmarked_launches"], 3)
            self.assertEqual(row["oom_launches"], 1)
            self.assertAlmostEqual(row["weighted_cpu_efficiency"], 0.75)
            connection.close()

    def test_requested_and_used_resources_are_kept_per_job(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            rows = connection.execute(
                "SELECT * FROM benchmark_job ORDER BY launch_index"
            ).fetchall()
            self.assertEqual(len(rows), 3)
            first = rows[0]
            self.assertEqual(first["run_id"], "20260825-101500")
            self.assertEqual(first["rule"], "assembly")
            self.assertEqual(first["requested_mem_mb"], 65536.0)
            self.assertEqual(first["max_rss_mb"], 65000.0)
            self.assertEqual(first["requested_runtime_min"], 720.0)
            self.assertEqual(first["elapsed_sec"], 1800.0)
            self.assertEqual(first["state"], "OUT_OF_MEMORY")
            # The TSV writes booleans as Python literals.
            self.assertEqual((first["oom"], first["timeout"]), (1, 0))
            self.assertEqual(rows[1]["attempt"], 2)
            connection.close()

    def test_rule_level_benchmark_is_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            rows = {
                row["rule"]: row
                for row in connection.execute("SELECT * FROM benchmark_rule")
            }
            self.assertEqual(set(rows), {"assembly", "binning"})
            self.assertEqual(rows["assembly"]["launches"], 2)
            self.assertEqual(rows["assembly"]["retries"], 1)
            self.assertEqual(rows["assembly"]["median_requested_mem_mb"], 98304.0)
            self.assertEqual(rows["binning"]["failed_launches"], 0)
            connection.close()

    def test_failures_are_recorded_once_per_rule_and_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root, sections=("resources",))
            rows = {
                (row["rule"], row["target"]): row
                for row in connection.execute("SELECT * FROM run_failure")
            }
            self.assertEqual(len(rows), 3)
            recovered = rows[("assembly", "A1")]
            self.assertEqual(recovered["run_id"], "20260825-101500")
            self.assertEqual(recovered["attempts"], 2)
            self.assertEqual(recovered["status"], "recovered")
            self.assertEqual(recovered["category"], "out-of-memory")
            self.assertIn("--memory-multiplier", recovered["action"])
            self.assertEqual(rows[("singlem", "S2")]["status"], "failed")
            # A workflow-level error belongs to no job, and keeps its own row.
            self.assertEqual(rows[("bowtie2", "workflow")]["category"], "missing-input")
            connection.close()

    def test_a_run_that_lost_no_job_leaves_the_failure_table_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            (root / "logging" / "drakkar_20260825-101500.failures.tsv").unlink()
            connection = self.build(root, sections=("resources",))
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM run_failure").fetchone()[0], 0
            )
            connection.close()

    def test_a_run_without_benchmark_files_still_ingests(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root / "drakkar_20260825-101500.yaml", """
                run_id: '20260825-101500'
                command: preprocessing
                status: completed
                """)
            connection = self.build(root, sections=("resources",))
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM run").fetchone()[0], 1
            )
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM benchmark_job").fetchone()[0],
                0,
            )
            connection.close()

    def test_ingest_log_records_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                'SELECT section, source_file, rows_ingested FROM ingest_log '
                'WHERE table_name = "genome_taxonomy"'
            ).fetchone()
            self.assertEqual(row["section"], "taxonomy")
            self.assertTrue(row["source_file"].endswith("genome_taxonomy.tsv"))
            self.assertEqual(row["rows_ingested"], 2)
            connection.close()

    def test_selected_sections_leave_other_tables_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root, sections=parse_sections("taxonomy"))
            taxonomy = connection.execute(
                "SELECT COUNT(*) FROM genome_taxonomy"
            ).fetchone()[0]
            samples = connection.execute("SELECT COUNT(*) FROM sample").fetchone()[0]
            self.assertEqual(taxonomy, 2)
            self.assertEqual(samples, 0)
            connection.close()

    def test_chunked_ingest_matches_single_pass_ingest(self):
        original = report_ingest.CHUNK_SIZE
        report_ingest.CHUNK_SIZE = 2
        try:
            with tempfile.TemporaryDirectory() as tmp:
                root = self.build_output_dir(Path(tmp))
                connection = self.build(root)
                annotations = connection.execute(
                    "SELECT COUNT(*) FROM gene_annotation"
                ).fetchone()[0]
                genes = connection.execute("SELECT COUNT(*) FROM gene").fetchone()[0]
                terms = connection.execute(
                    "SELECT COUNT(*) FROM annotation_term"
                ).fetchone()[0]
                self.assertEqual(annotations, 3)
                self.assertEqual(genes, 3)
                self.assertEqual(terms, 2)
                connection.close()
        finally:
            report_ingest.CHUNK_SIZE = original


class LegacyLayoutIngestTests(ReportFixtureMixin, unittest.TestCase):
    """Output directories written before the ``logging/`` directory existed.

    Their run metadata sits in the output root, their benchmark roll-up beside
    it and their benchmark tables in ``benchmark/``. The report still finds all
    of it, so re-reporting an older directory does not lose its runs.
    """

    def build(self, root, sections=None):
        db_path = Path(root) / "drakkar.db"
        report_command.build_database(root, sections or ("resources",), db_path)
        return connect(db_path)

    def test_run_metadata_in_the_output_root_is_ingested(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_legacy_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute("SELECT * FROM run").fetchone()
            self.assertEqual(row["run_id"], "20260825-101500")
            self.assertEqual(row["status"], "completed")
            self.assertEqual(row["drakkar_version"], "2.0.0")
            connection.close()

    def test_legacy_benchmark_artefacts_are_ingested(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_legacy_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute("SELECT * FROM run_benchmark").fetchone()
            self.assertEqual(row["run_id"], "20260825-101500")
            self.assertEqual(row["benchmarked_launches"], 3)
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM benchmark_job").fetchone()[0], 3
            )
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM benchmark_rule").fetchone()[0], 2
            )
            connection.close()

    def test_benchmark_summary_does_not_overwrite_the_run_row(self):
        # In this layout drakkar_<run_id>_resources.yaml sits in the output root
        # and matches the run metadata glob, carrying the same run_id, so it
        # must not be read as a run itself. The current layout nests it under
        # logging/benchmark/, where the glob cannot reach it.
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_legacy_output_dir(Path(tmp))
            connection = self.build(root)
            rows = connection.execute("SELECT * FROM run").fetchall()
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["status"], "completed")
            self.assertEqual(rows[0]["drakkar_version"], "2.0.0")
            connection.close()

    def test_legacy_failure_table_is_ingested(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_legacy_output_dir(Path(tmp))
            connection = self.build(root)
            rows = connection.execute(
                "SELECT run_id, rule, attempts FROM run_failure ORDER BY rule"
            ).fetchall()
            self.assertEqual([row["rule"] for row in rows],
                             ["assembly", "bowtie2", "singlem"])
            self.assertEqual(rows[0]["run_id"], "20260825-101500")
            self.assertEqual(rows[0]["attempts"], 2)
            connection.close()

    def test_probe_finds_the_legacy_resources_section(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_legacy_output_dir(Path(tmp))
            entry = {item["section"]: item for item in probe(root)}["resources"]
            self.assertTrue(entry["available"])
            self.assertEqual(entry["missing"], [])
            self.assertIn("drakkar_20260825-101500.yaml", entry["present"])
            self.assertIn("drakkar_20260825-101500_resources.yaml", entry["present"])


class AmrIngestTests(ReportFixtureMixin, unittest.TestCase):
    """The assembly-level AMR tables written by `aggregate_amr`."""

    def build(self, root):
        db_path = Path(root) / "drakkar.db"
        report_command.build_database(root, parse_sections(None), db_path)
        connection = sqlite3.connect(db_path)
        connection.row_factory = sqlite3.Row
        self.addCleanup(connection.close)
        return connection

    def test_the_summary_and_the_qc_roll_up_land_in_one_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            row = self.build(root).execute(
                "SELECT * FROM amr_assembly WHERE assembly_id = 'A1'"
            ).fetchone()
            # From assembly_summary.tsv.
            self.assertEqual(row["assembly_type"], "metagenome")
            self.assertEqual(row["contig_count"], 1200)
            self.assertEqual(row["amr_loci"], 3)
            # From amr_qc.tsv, which the summary does not carry.
            self.assertEqual(row["rgi_hits_without_coordinates"], 1)
            self.assertEqual(row["mobility_links"], 2)

    def test_loci_keep_the_caller_support_and_the_agreement_between_them(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            rows = self.build(root).execute(
                "SELECT locus_id, support_status, concordance, source_count "
                "FROM amr_locus WHERE assembly_id = 'A1' ORDER BY locus_id"
            ).fetchall()
            self.assertEqual(
                [(row["locus_id"], row["support_status"], row["concordance"])
                 for row in rows],
                [("L1", "amrfinder_and_rgi", "gene_match"),
                 ("L2", "amrfinder_only", "single_source"),
                 ("L3", "rgi_only", "single_source")],
            )
            self.assertEqual(rows[0]["source_count"], 2)

    def test_hits_keep_their_caller_and_drop_the_raw_details_blob(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            row = connection.execute(
                "SELECT source, identity, detection_grade, is_partial "
                "FROM amr_hit WHERE hit_id = 'H2'"
            ).fetchone()
            self.assertEqual(row["source"], "rgi")
            self.assertAlmostEqual(row["identity"], 99.6)
            self.assertEqual(row["detection_grade"], "Perfect")
            self.assertEqual(row["is_partial"], 0)
            columns = {
                item[1] for item in
                connection.execute("PRAGMA table_info(amr_hit)").fetchall()
            }
            self.assertNotIn("raw_details", columns)

    def test_a_partial_hit_is_stored_as_a_flag_not_as_its_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            row = self.build(root).execute(
                "SELECT is_partial FROM amr_hit WHERE hit_id = 'H3'"
            ).fetchone()
            self.assertEqual(row["is_partial"], 1)

    def test_drug_classes_carry_one_row_per_hit_and_class(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            rows = self.build(root).execute(
                "SELECT drug_class, COUNT(*) AS hits FROM amr_drug_class "
                "GROUP BY drug_class ORDER BY drug_class"
            ).fetchall()
            self.assertEqual(
                [(row["drug_class"], row["hits"]) for row in rows],
                [("beta-lactam", 3), ("quinolone", 1), ("sulfonamide", 1)],
            )

    def test_mobility_links_the_loci_to_the_regions_they_sit_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            rows = self.build(root).execute(
                "SELECT m.locus_id, m.context_type, r.topology, r.score "
                "FROM amr_mobility AS m "
                "JOIN amr_mobility_region AS r "
                "  ON r.assembly_id = m.assembly_id AND r.region_id = m.region_id "
                "ORDER BY m.locus_id"
            ).fetchall()
            self.assertEqual(
                [(row["locus_id"], row["context_type"], row["topology"])
                 for row in rows],
                [("L1", "plasmid", "circular"), ("L3", "provirus", "linear")],
            )
            self.assertAlmostEqual(rows[0]["score"], 0.94)

    def test_re_ingesting_the_same_directory_does_not_duplicate_rows(self):
        # `drakkar reporting` reuses a matching database rather than rebuilding
        # it, so every AMR loader has to be safe to run over its own output.
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            db_path = root / "drakkar.db"
            report_command.build_database(root, parse_sections(None), db_path)
            connection = sqlite3.connect(db_path)
            self.addCleanup(connection.close)
            tables = ("amr_assembly", "amr_hit", "amr_locus", "amr_drug_class",
                      "amr_mobility_region", "amr_mobility")
            first = {
                table: connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                for table in tables
            }
            report_command.build_database(root, parse_sections(None), db_path)
            for table in tables:
                with self.subTest(table=table):
                    self.assertEqual(
                        connection.execute(
                            f"SELECT COUNT(*) FROM {table}"
                        ).fetchone()[0],
                        first[table],
                    )

    def test_an_output_directory_without_amr_leaves_the_tables_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root / "preprocessing.tsv",
                  "sample\treads_pre_fastp\nS1\t1000\n")
            connection = self.build(root)
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM amr_assembly").fetchone()[0], 0
            )
            self.assertEqual(
                connection.execute(
                    "SELECT COUNT(*) FROM ingest_log WHERE section = 'amr'"
                ).fetchone()[0],
                0,
            )


class BinOriginTests(ReportFixtureMixin, unittest.TestCase):
    """The per-assembly Binette reports behind `cataloging/final`."""

    def build(self, root):
        db_path = Path(root) / "drakkar.db"
        report_command.build_database(root, parse_sections(None), db_path)
        connection = sqlite3.connect(db_path)
        connection.row_factory = sqlite3.Row
        self.addCleanup(connection.close)
        return connection

    def test_bins_are_named_the_way_drakkar_renames_them(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            names = [
                row[0] for row in
                self.build(root).execute("SELECT bin_name FROM assembly_bin ORDER BY bin_name")
            ]
            self.assertEqual(
                names, ["A1_bin_1", "A1_bin_2", "A1_bin_3", "A2_bin_1"]
            )

    def test_binner_names_are_normalized_to_the_ones_the_counts_use(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            binners = [
                row[0] for row in connection.execute(
                    "SELECT DISTINCT binner FROM assembly_bin_origin ORDER BY binner"
                )
            ]
            # Binette names the sets after its input directories — `metabat`,
            # `semibin` — while the counts table uses the tool versions.
            self.assertEqual(binners, ["binette", "maxbin2", "metabat2", "semibin2"])

    def test_a_bin_several_binners_produced_names_all_of_them(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            binners = [
                row[0] for row in connection.execute(
                    "SELECT binner FROM assembly_bin_origin "
                    "WHERE bin_name = 'A1_bin_2' ORDER BY binner"
                )
            ]
            self.assertEqual(binners, ["metabat2", "semibin2"])

    def test_a_bin_binette_built_is_not_original(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            row = self.build(root).execute(
                "SELECT is_original, origin FROM assembly_bin WHERE bin_name = 'A1_bin_3'"
            ).fetchone()
            self.assertEqual(row["is_original"], 0)
            self.assertEqual(row["origin"], "binette")

    def test_absent_reports_leave_the_rest_of_cataloging_intact(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root / "cataloging.tsv", """
                assembly\tmetabat2_bins\tfinal_bins
                A1\t5\t3
                """)
            connection = self.build(root)
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM assembly").fetchone()[0], 1
            )
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM assembly_bin").fetchone()[0], 0
            )


class DereplicationClusterTests(ReportFixtureMixin, unittest.TestCase):
    """dRep's cluster assignments and pairwise comparisons."""

    def build(self, root):
        db_path = Path(root) / "drakkar.db"
        report_command.build_database(root, parse_sections(None), db_path)
        connection = sqlite3.connect(db_path)
        connection.row_factory = sqlite3.Row
        self.addCleanup(connection.close)
        return connection

    def test_cluster_assignments_and_winners_are_loaded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            rows = {
                row["genome"]: row for row in
                connection.execute("SELECT * FROM genome_cluster")
            }
            self.assertEqual(len(rows), 5)
            self.assertEqual(rows["A1_bin_1"]["secondary_cluster"], "1_1")
            self.assertEqual(rows["A1_bin_1"]["is_representative"], 1)
            # Clustered with A1_bin_1 and lost to it.
            self.assertEqual(rows["A1_bin_2"]["is_representative"], 0)

    def test_only_the_nearest_mash_neighbour_is_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            connection = self.build(root)
            rows = {
                row["genome"]: row["nearest_mash_distance"] for row in
                connection.execute(
                    "SELECT genome, nearest_mash_distance FROM genome_cluster"
                )
            }
            # A1_bin_1 is 0.02 from A1_bin_2 and 0.41 from A1_bin_3; the
            # self-comparison at distance 0 must not win.
            self.assertAlmostEqual(rows["A1_bin_1"], 0.02)
            self.assertAlmostEqual(rows["A2_bin_2"], 0.45)

    def test_a_bin_with_no_near_neighbour_is_distinguishable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            near = self.build(root).execute(
                "SELECT COUNT(*) FROM genome_cluster WHERE nearest_mash_distance <= 0.1"
            ).fetchone()[0]
            self.assertEqual(near, 4)

    def test_reciprocal_comparisons_collapse_to_one_averaged_pair(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            rows = self.build(root).execute(
                "SELECT genome_a, genome_b, ani, primary_cluster "
                "FROM genome_comparison ORDER BY genome_a"
            ).fetchall()
            self.assertEqual(len(rows), 2)
            self.assertEqual((rows[0]["genome_a"], rows[0]["genome_b"]),
                             ("A1_bin_1", "A1_bin_2"))
            # dRep clusters on the average of the two directions.
            self.assertAlmostEqual(rows[0]["ani"], 0.998)
            self.assertEqual(rows[0]["primary_cluster"], "1")

    def test_absent_drep_tables_leave_the_summary_intact(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root / "dereplicating.tsv", """
                input_bin_number\tdereplication_ani\toutput_bin_number
                11\t0.98\t6
                """)
            connection = self.build(root)
            self.assertEqual(
                connection.execute(
                    "SELECT input_bin_number FROM dereplication"
                ).fetchone()[0], 11
            )
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM genome_cluster").fetchone()[0], 0
            )


class TaxonomyTreeTests(unittest.TestCase):
    """The pruned GTDB-Tk placement trees, stored beside the lineages."""

    SUMMARY = """
        user_genome\tclassification\tclassification_method
        MAG_1\td__Bacteria;p__Bacillota;c__Clostridia\ttopology and ANI
        MAG_2\td__Bacteria;p__Bacteroidota;c__Bacteroidia\ttopology and ANI
        MAG_3\td__Archaea;p__Methanobacteriota;c__Methanobacteria\ttopology and ANI
        MAG_4\td__Archaea;p__Methanobacteriota;c__Methanobacteria\ttopology and ANI
        """

    def ingest(self, root, bacteria=None, archaea=None):
        write(root / "annotating/genome_taxonomy.tsv", self.SUMMARY)
        if bacteria is not None:
            write(root / "annotating/bacteria.tree", bacteria)
        if archaea is not None:
            write(root / "annotating/archaea.tree", archaea)
        connection = connect(root / "drakkar.db")
        create_schema(connection)
        report_ingest.ingest_taxonomy(connection, root)
        return connection

    def trees(self, connection):
        return {
            row[0]: (row[1], row[2])
            for row in connection.execute(
                "SELECT domain, newick, tip_count FROM genome_tree"
            )
        }

    def test_each_domain_tree_is_stored_with_its_tip_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            connection = self.ingest(
                root,
                bacteria="(MAG_1:0.12,MAG_2:0.09);\n",
                archaea="(MAG_3:0.20,MAG_4:0.18);\n",
            )
            trees = self.trees(connection)
            self.assertEqual(set(trees), {"bacteria", "archaea"})
            self.assertEqual(trees["bacteria"][0], "(MAG_1:0.12,MAG_2:0.09);")
            self.assertEqual(trees["bacteria"][1], 2)
            self.assertEqual(trees["archaea"][1], 2)
            connection.close()

    def test_a_catalogue_without_archaea_stores_only_the_bacterial_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            connection = self.ingest(root, bacteria="(MAG_1:0.12,MAG_2:0.09);")
            self.assertEqual(set(self.trees(connection)), {"bacteria"})
            connection.close()

    def test_the_empty_tree_the_pruning_rule_writes_is_not_stored(self):
        # `gtdbtk_pruned_trees` touches an empty bacteria.tree when GTDB-Tk
        # produced no classify tree at all.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "annotating").mkdir(parents=True)
            (root / "annotating/bacteria.tree").write_text("", encoding="utf-8")
            connection = self.ingest(root)
            self.assertEqual(self.trees(connection), {})
            connection.close()

    def test_a_tree_of_one_tip_has_no_topology_and_is_not_stored(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            connection = self.ingest(root, bacteria="MAG_1:0.1;")
            self.assertEqual(self.trees(connection), {})
            connection.close()

    def test_an_unparsable_tree_is_skipped_rather_than_raised(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            connection = self.ingest(root, bacteria="((MAG_1,MAG_2);")
            self.assertEqual(self.trees(connection), {})
            connection.close()

    def test_the_taxonomy_lineages_are_ingested_either_way(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            connection = self.ingest(root)
            self.assertEqual(
                connection.execute(
                    "SELECT COUNT(*) FROM genome_taxonomy"
                ).fetchone()[0],
                4,
            )
            connection.close()

    def test_each_tree_is_named_in_the_ingest_log_with_its_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            connection = self.ingest(
                root, bacteria="(MAG_1:0.1,MAG_2:0.1);",
                archaea="(MAG_3:0.1,MAG_4:0.1);",
            )
            rows = {
                row[0]: row[1]
                for row in connection.execute(
                    "SELECT table_name, source_file FROM ingest_log "
                    "WHERE table_name LIKE 'genome_tree%'"
                )
            }
            self.assertEqual(
                set(rows), {"genome_tree:bacteria", "genome_tree:archaea"}
            )
            self.assertTrue(rows["genome_tree:bacteria"].endswith("bacteria.tree"))
            connection.close()

    def test_the_probe_names_the_trees_it_found(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root / "annotating/genome_taxonomy.tsv", self.SUMMARY)
            write(root / "annotating/bacteria.tree", "(MAG_1:0.1,MAG_2:0.1);")
            entry = {item["section"]: item for item in probe(root)}["taxonomy"]
            self.assertTrue(entry["available"])
            self.assertIn("annotating/bacteria.tree", entry["present"])
            self.assertNotIn("annotating/archaea.tree", entry["present"])

    def test_a_taxonomy_without_trees_is_still_available(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write(root / "annotating/genome_taxonomy.tsv", self.SUMMARY)
            entry = {item["section"]: item for item in probe(root)}["taxonomy"]
            self.assertTrue(entry["available"])
            self.assertEqual(entry["missing"], [])


class MissingSourceTests(unittest.TestCase):
    def test_loaders_return_none_when_sources_are_absent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            connection = connect(root / "drakkar.db")
            create_schema(connection)
            for name, loader in report_ingest.SECTION_LOADERS.items():
                self.assertIsNone(loader(connection, root), name)
            connection.close()


class ReportCommandTests(ReportFixtureMixin, unittest.TestCase):
    def test_missing_output_directory_fails(self):
        self.assertEqual(report_command.run_report("/nonexistent-drakkar-dir"), 1)

    def test_unknown_section_fails_without_building(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            self.assertEqual(report_command.run_report(root, sections="bogus"), 1)
            self.assertFalse(report_command.database_path(root).exists())

    def test_empty_directory_fails_with_no_sections(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(report_command.run_report(tmp, db_only=True), 1)

    def test_successful_build_writes_the_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            self.assertEqual(report_command.run_report(root, db_only=True), 0)
            self.assertTrue(report_command.database_path(root).exists())

    def test_schema_mismatch_blocks_until_forced(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            self.assertEqual(report_command.run_report(root, db_only=True), 0)
            connection = sqlite3.connect(report_command.database_path(root))
            connection.execute("UPDATE schema_version SET version = 99")
            connection.commit()
            connection.close()
            self.assertEqual(report_command.run_report(root, db_only=True), 1)
            self.assertEqual(
                report_command.run_report(root, db_only=True, force=True), 0
            )

    def test_a_database_at_the_output_root_is_moved_into_reporting(self):
        # Directories built before the reporting directory existed keep their
        # database at the root; it is moved rather than rebuilt beside itself.
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            self.assertEqual(report_command.run_report(root, db_only=True), 0)
            legacy = report_command.legacy_database_path(root)
            report_command.database_path(root).replace(legacy)
            (root / report_command.REPORTING_DIRNAME).rmdir()

            self.assertEqual(report_command.run_report(root, db_only=True), 0)
            self.assertFalse(legacy.exists())
            self.assertTrue(report_command.database_path(root).exists())

    def test_a_report_at_the_output_root_is_still_found(self):
        # Moving where reports are written must not hide the ones already
        # rendered next to the source tables.
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            older = root / report_command.report_name()
            older.write_text("<html></html>", encoding="utf-8")
            self.assertIn(older, report_command.find_reports(root))

    def test_rebuild_is_idempotent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build_output_dir(Path(tmp))
            report_command.run_report(root, db_only=True)
            report_command.run_report(root, db_only=True, force=True)
            connection = sqlite3.connect(report_command.database_path(root))
            count = connection.execute("SELECT COUNT(*) FROM gene_annotation").fetchone()[0]
            self.assertEqual(count, 3)
            connection.close()


if __name__ == "__main__":
    unittest.main()
