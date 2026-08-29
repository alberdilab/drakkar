"""Assembly-level AMRFinderPlus, CARD/RGI, and geNomad workflow."""

import shlex
from pathlib import Path


PACKAGE_DIR = config["package_dir"]
AMRFINDER_DB = config["AMRFINDER_DB"]
CARD_DB = config["CARD_DB"]
GENOMAD_DB = config["GENOMAD_DB"]
GENOMAD_MODULE = config["GENOMAD_MODULE"]
GENOMAD_SPLITS = positive_config_int("genomad_splits", config.get("GENOMAD_SPLITS", 8))
AMR_LOCUS_OVERLAP = float(config.get("locus_overlap", 0.8))
RGI_ALIGNMENT_TOOL = str(config.get("rgi_alignment_tool", "DIAMOND")).upper()
GENOMAD_PRESET = str(config.get("genomad_preset", "default")).lower()
DRAKKAR_VERSION = str(config.get("drakkar_version", "unknown"))


def amr_config_bool(key, default=False):
    value = config.get(key, default)
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


RGI_INCLUDE_LOOSE = amr_config_bool("rgi_include_loose")
RGI_INCLUDE_NUDGE = amr_config_bool("rgi_include_nudge")


def assembly_metadata(wildcards):
    return AMR_ASSEMBLIES[wildcards.assembly]


def prodigal_mode(wildcards):
    return "single" if assembly_metadata(wildcards)["assembly_type"] == "isolate" else "meta"


def amrfinder_organism_flag(wildcards):
    organism = str(assembly_metadata(wildcards).get("organism") or "").strip()
    return f"--organism {shlex.quote(organism)}" if organism else ""


def rgi_low_quality_flag(wildcards):
    return "--low_quality" if assembly_metadata(wildcards)["assembly_type"] == "metagenome" else ""


def genomad_preset_flag():
    return "" if GENOMAD_PRESET == "default" else f"--{GENOMAD_PRESET}"


rule stage_amr_assembly:
    input:
        lambda wildcards: assembly_metadata(wildcards)["path"]
    output:
        f"{OUTPUT_DIR}/amr/staged/{{assembly}}.fna"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(1024 * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(10, int(input.size_mb * 2)) * 2 ** (attempt - 1))
    shell:
        r"""
        mkdir -p $(dirname {output:q})
        case {input:q} in
            *.gz) gzip -dc {input:q} > {output:q} ;;
            *) cp {input:q} {output:q} ;;
        esac
        """


rule amr_prodigal:
    input:
        assembly=f"{OUTPUT_DIR}/amr/staged/{{assembly}}.fna"
    output:
        proteins=f"{OUTPUT_DIR}/amr/raw/prodigal/{{assembly}}.faa",
        genes=f"{OUTPUT_DIR}/amr/raw/prodigal/{{assembly}}.ffn",
        gff=f"{OUTPUT_DIR}/amr/raw/prodigal/{{assembly}}.gff",
        amrfinder_gff=f"{OUTPUT_DIR}/amr/raw/prodigal/{{assembly}}.amrfinder.gff"
    params:
        mode=prodigal_mode,
        converter=f"{PACKAGE_DIR}/workflow/scripts/prepare_amrfinder_gff.py"
    threads: 1
    conda:
        f"{PACKAGE_DIR}/workflow/envs/amr_amrfinder.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(2048, int(input.size_mb * 8)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb * 3)) * 2 ** (attempt - 1))
    shell:
        r"""
        mkdir -p $(dirname {output.proteins:q})
        prodigal \
            -i {input.assembly:q} \
            -a {output.proteins:q} \
            -d {output.genes:q} \
            -f gff \
            -o {output.gff:q} \
            -p {params.mode:q}
        python {params.converter:q} {output.gff:q} {output.amrfinder_gff:q} --proteins {output.proteins:q}
        """


rule amrfinderplus:
    input:
        assembly=f"{OUTPUT_DIR}/amr/staged/{{assembly}}.fna",
        proteins=f"{OUTPUT_DIR}/amr/raw/prodigal/{{assembly}}.faa",
        gff=f"{OUTPUT_DIR}/amr/raw/prodigal/{{assembly}}.amrfinder.gff"
    output:
        f"{OUTPUT_DIR}/amr/raw/amrfinder/{{assembly}}.tsv"
    params:
        database=AMRFINDER_DB,
        organism=amrfinder_organism_flag
    threads: 4
    conda:
        f"{PACKAGE_DIR}/workflow/envs/amr_amrfinder.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(4096, int(input.size_mb * 20)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(30, int(input.size_mb * 5)) * 2 ** (attempt - 1))
    shell:
        r"""
        mkdir -p $(dirname {output:q})
        amrfinder \
            --nucleotide {input.assembly:q} \
            --protein {input.proteins:q} \
            --gff {input.gff:q} \
            --annotation_format standard \
            --database {params.database:q} \
            --threads {threads} \
            --print_node \
            --name {wildcards.assembly:q} \
            {params.organism} \
            --output {output:q}
        """


rule rgi_card:
    input:
        assembly=f"{OUTPUT_DIR}/amr/staged/{{assembly}}.fna"
    output:
        table=f"{OUTPUT_DIR}/amr/raw/rgi/{{assembly}}.txt",
        json=f"{OUTPUT_DIR}/amr/raw/rgi/{{assembly}}.json"
    params:
        card_root=CARD_DB,
        output_base=lambda wildcards: str(
            Path(f"{OUTPUT_DIR}/amr/raw/rgi/{wildcards.assembly}").resolve()
        ),
        alignment=RGI_ALIGNMENT_TOOL,
        low_quality=rgi_low_quality_flag,
        loose="--include_loose" if RGI_INCLUDE_LOOSE else "",
        nudge="--include_nudge" if RGI_INCLUDE_NUDGE else ""
    threads: 8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/amr_rgi.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8192, int(input.size_mb * 40)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(45, int(input.size_mb * 10)) * 2 ** (attempt - 1))
    shell:
        r"""
        mkdir -p $(dirname {output.table:q})
        cd {params.card_root:q}
        rgi main \
            --input_sequence {input.assembly:q} \
            --output_file {params.output_base:q} \
            --input_type contig \
            --alignment_tool {params.alignment:q} \
            --threads {threads} \
            --local \
            --clean \
            {params.low_quality} \
            {params.loose} \
            {params.nudge}
        test -s {output.table:q}
        test -s {output.json:q}
        """


rule genomad_amr_context:
    input:
        assembly=f"{OUTPUT_DIR}/amr/staged/{{assembly}}.fna"
    output:
        plasmid_summary=f"{OUTPUT_DIR}/amr/raw/genomad/{{assembly}}/{{assembly}}_summary/{{assembly}}_plasmid_summary.tsv",
        plasmid_genes=f"{OUTPUT_DIR}/amr/raw/genomad/{{assembly}}/{{assembly}}_summary/{{assembly}}_plasmid_genes.tsv",
        virus_summary=f"{OUTPUT_DIR}/amr/raw/genomad/{{assembly}}/{{assembly}}_summary/{{assembly}}_virus_summary.tsv",
        virus_genes=f"{OUTPUT_DIR}/amr/raw/genomad/{{assembly}}/{{assembly}}_summary/{{assembly}}_virus_genes.tsv"
    params:
        module=GENOMAD_MODULE,
        outdir=f"{OUTPUT_DIR}/amr/raw/genomad/{{assembly}}",
        database=GENOMAD_DB,
        preset=genomad_preset_flag()
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(16 * 1024, int(input.size_mb * 1024 * 50)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(30, int(input.size_mb * 20)) * 2 ** (attempt - 1))
    shell:
        r"""
        module purge
        module load {params.module}
        unset PYTHONPATH
        export PYTHONNOUSERSITE=1
        genomad end-to-end \
            --cleanup \
            --splits {GENOMAD_SPLITS} \
            --threads {threads} \
            {params.preset} \
            {input.assembly:q} \
            {params.outdir:q} \
            {params.database:q}
        """


rule digest_amr_assembly:
    input:
        amrfinder=f"{OUTPUT_DIR}/amr/raw/amrfinder/{{assembly}}.tsv",
        rgi=f"{OUTPUT_DIR}/amr/raw/rgi/{{assembly}}.txt",
        plasmid_summary=f"{OUTPUT_DIR}/amr/raw/genomad/{{assembly}}/{{assembly}}_summary/{{assembly}}_plasmid_summary.tsv",
        plasmid_genes=f"{OUTPUT_DIR}/amr/raw/genomad/{{assembly}}/{{assembly}}_summary/{{assembly}}_plasmid_genes.tsv",
        virus_summary=f"{OUTPUT_DIR}/amr/raw/genomad/{{assembly}}/{{assembly}}_summary/{{assembly}}_virus_summary.tsv",
        virus_genes=f"{OUTPUT_DIR}/amr/raw/genomad/{{assembly}}/{{assembly}}_summary/{{assembly}}_virus_genes.tsv"
    output:
        hits=f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/amr_hits.tsv",
        loci=f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/amr_loci.tsv",
        drugs=f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/amr_drug_classes.tsv",
        regions=f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/mobility_regions.tsv",
        mobility=f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/amr_mobility.tsv",
        digest=f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/digest.tsv",
        qc=f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/qc.json"
    params:
        script=f"{PACKAGE_DIR}/workflow/scripts/amr_digest.py"
    threads: 1
    conda:
        f"{PACKAGE_DIR}/workflow/envs/amr_digest.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(4096 * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(60 * 2 ** (attempt - 1))
    shell:
        r"""
        python {params.script:q} \
            --assembly-id {wildcards.assembly:q} \
            --amrfinder {input.amrfinder:q} \
            --rgi {input.rgi:q} \
            --plasmid-summary {input.plasmid_summary:q} \
            --plasmid-genes {input.plasmid_genes:q} \
            --virus-summary {input.virus_summary:q} \
            --virus-genes {input.virus_genes:q} \
            --minimum-overlap {AMR_LOCUS_OVERLAP} \
            --hits-output {output.hits:q} \
            --loci-output {output.loci:q} \
            --drugs-output {output.drugs:q} \
            --regions-output {output.regions:q} \
            --mobility-output {output.mobility:q} \
            --digest-output {output.digest:q} \
            --qc-output {output.qc:q}
        """


rule aggregate_amr:
    input:
        manifest=f"{OUTPUT_DIR}/data/amr_assemblies.json",
        hits=expand(f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/amr_hits.tsv", assembly=assemblies),
        loci=expand(f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/amr_loci.tsv", assembly=assemblies),
        drugs=expand(f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/amr_drug_classes.tsv", assembly=assemblies),
        regions=expand(f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/mobility_regions.tsv", assembly=assemblies),
        mobility=expand(f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/amr_mobility.tsv", assembly=assemblies),
        digests=expand(f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/digest.tsv", assembly=assemblies),
        qc=expand(f"{OUTPUT_DIR}/amr/assemblies/{{assembly}}/qc.json", assembly=assemblies)
    output:
        hits=f"{OUTPUT_DIR}/amr/amr_hits.tsv.xz",
        loci=f"{OUTPUT_DIR}/amr/amr_loci.tsv.xz",
        drugs=f"{OUTPUT_DIR}/amr/amr_drug_classes.tsv.xz",
        regions=f"{OUTPUT_DIR}/amr/mobility_regions.tsv.xz",
        mobility=f"{OUTPUT_DIR}/amr/amr_mobility.tsv.xz",
        summary=f"{OUTPUT_DIR}/amr/assembly_summary.tsv",
        qc=f"{OUTPUT_DIR}/amr/amr_qc.tsv",
        manifest=f"{OUTPUT_DIR}/amr/manifest.yaml"
    params:
        script=f"{PACKAGE_DIR}/workflow/scripts/aggregate_amr.py",
        assembly_dir=f"{OUTPUT_DIR}/amr/assemblies",
        output_dir=f"{OUTPUT_DIR}/amr",
        loose="--rgi-include-loose" if RGI_INCLUDE_LOOSE else "",
        nudge="--rgi-include-nudge" if RGI_INCLUDE_NUDGE else "",
        genomad_version=str(GENOMAD_MODULE).rsplit("/", 1)[-1]
    threads: 1
    conda:
        f"{PACKAGE_DIR}/workflow/envs/amr_digest.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(4096, int(input.size_mb * 4)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(30, len(assemblies) * 2) * 2 ** (attempt - 1))
    shell:
        r"""
        python {params.script:q} \
            --input-manifest {input.manifest:q} \
            --per-assembly-dir {params.assembly_dir:q} \
            --output-dir {params.output_dir:q} \
            --minimum-overlap {AMR_LOCUS_OVERLAP} \
            --rgi-alignment-tool {RGI_ALIGNMENT_TOOL:q} \
            {params.loose} \
            {params.nudge} \
            --genomad-preset {GENOMAD_PRESET:q} \
            --genomad-splits {GENOMAD_SPLITS} \
            --amrfinder-db {AMRFINDER_DB:q} \
            --card-db {CARD_DB:q} \
            --genomad-db {GENOMAD_DB:q} \
            --drakkar-version {DRAKKAR_VERSION:q} \
            --genomad-version {params.genomad_version:q}
        """
