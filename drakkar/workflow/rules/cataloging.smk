####
# Config variables
####

PACKAGE_DIR = config["package_dir"]

# Software modules
MEGAHIT_MODULE = config["MEGAHIT_MODULE"]
BOWTIE2_MODULE = config["BOWTIE2_MODULE"]
SAMTOOLS_MODULE = config["SAMTOOLS_MODULE"]
QUAST_MODULE = config["QUAST_MODULE"]
METABAT2_MODULE = config["METABAT2_MODULE"]
MAXBIN2_MODULE = config["MAXBIN2_MODULE"]
FRAGGENESCAN_MODULE = config["FRAGGENESCAN_MODULE"]
HMMER_MODULE = config["HMMER_MODULE"]
DIAMOND_MODULE = config["DIAMOND_MODULE"]
CHECKM2_MODULE = config["CHECKM2_MODULE"]
BINETTE_MODULE = config["BINETTE_MODULE"]
MIN_COMPLETENESS = int(config.get("MIN_COMPLETENESS", 70))
MAX_CONTAMINATION = int(config.get("MAX_CONTAMINATION", 10))
MIN_BIN_LENGTH = int(config.get("MIN_BIN_LENGTH", 200000))
MAX_BIN_LENGTH = int(config.get("MAX_BIN_LENGTH", 10000000))
# Assemblies below this size (MB of FASTA) are not worth binning with the
# heavier binners, and COMEBin cannot bin them at all: its marker-gene stage
# aborts when FragGeneScan/HMMsearch find no single-copy markers to seed the
# clustering. Such assemblies export an empty contig-to-bin table instead.
MIN_BINNING_ASSEMBLY_MB = int(config.get("MIN_BINNING_ASSEMBLY_MB", 10))

# Binette >=1.2 names the selected bins <prefix>_bin<n> and writes both the
# quality report and the bin FASTA files using that name.
BINETTE_REPORT_HEADER = "\\t".join(BINETTE_REPORT_COLUMNS) + "\\n"

# Databases
CHECKM2_DB = config["CHECKM2_DB"]

BINNER_ORDER = ("metabat", "maxbin", "semibin", "comebin")
BINNER_ALIASES = {
    "metabat": "metabat",
    "metabat2": "metabat",
    "maxbin": "maxbin",
    "maxbin2": "maxbin",
    "semibin": "semibin",
    "semibin2": "semibin",
    "comebin": "comebin",
}
BINNER_RULE_NAMES = {
    "metabat": "metabat2",
    "maxbin": "maxbin2",
    "semibin": "semibin2",
    "comebin": "comebin",
}


def normalize_binner_config(raw_binners):
    if raw_binners is None:
        return list(BINNER_ORDER)
    if isinstance(raw_binners, str):
        items = [item.strip().lower() for item in raw_binners.split(",") if item.strip()]
    else:
        items = [str(item).strip().lower() for item in raw_binners if str(item).strip()]
    if not items or "all" in items:
        return list(BINNER_ORDER)

    invalid = [item for item in items if item not in BINNER_ALIASES]
    if invalid:
        raise ValueError(
            f"Unsupported binner(s): {', '.join(invalid)}. "
            f"Options are: {', '.join(BINNER_ORDER)}."
        )

    selected = {BINNER_ALIASES[item] for item in items}
    return [binner for binner in BINNER_ORDER if binner in selected]


SELECTED_BINNERS = normalize_binner_config(config.get("binners"))


def binner_table_path(assembly, binner):
    rule_name = BINNER_RULE_NAMES[binner]
    return f"{OUTPUT_DIR}/cataloging/{rule_name}/{assembly}/{assembly}.tsv"


def selected_binner_tables(wildcards):
    return [binner_table_path(wildcards.assembly, binner) for binner in SELECTED_BINNERS]


def selected_binner_expand(binner):
    if binner not in SELECTED_BINNERS:
        return []
    return expand(binner_table_path("{assembly}", binner), assembly=assemblies)

####
# Workflow rules
####

rule assembly:
    input:
        r1=lambda wildcards: [path for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly] for path in PREPROCESSED_TO_READS1[sample]],
        r2=lambda wildcards: [path for sample in ASSEMBLY_TO_SAMPLES[wildcards.assembly] for path in PREPROCESSED_TO_READS2[sample]]
    output:
        f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna"
    params:
        megahit_module={MEGAHIT_MODULE},
        outputdir=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(min(1020*1024,max(8*1024, int(input.size_mb * 10)) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 20)) * 2 ** (attempt - 1))
    message: "Assembling {wildcards.assembly}..."
    shell:
        """
        module purge
        module load {params.megahit_module}
        rm -rf {params.outputdir}

        # Convert input list to a comma-separated string
        R1_FILES=$(echo {input.r1} | tr ' ' ',')
        R2_FILES=$(echo {input.r2} | tr ' ' ',')

        megahit \
            -t {threads} \
            --verbose \
            --min-contig-len 1500 \
            -1 $R1_FILES -2 $R2_FILES \
            -o {params.outputdir}
        mv {params.outputdir}/final.contigs.fa {params.outputdir}/final.contigs.raw.fa
        awk -v a="{wildcards.assembly}" 'BEGIN{{i=1}} /^>/{{print ">" a "_" i++; next}} {{print}}' \
            {params.outputdir}/final.contigs.raw.fa > {output}
        """

rule assembly_index:
    input:
        f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna"
    output:
        index=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.rev.2.bt2"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        basename=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 10)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 5)) * 2 ** (attempt - 1))
    message: "Indexing assembly {wildcards.assembly}..."
    shell:
        """
        if [ ! -s {input} ]; then
            echo "Assembly is empty, skipping bowtie2-build..."
            touch {output.index}
        else
            module purge
            module load {params.bowtie2_module}
            bowtie2-build {input} {params.basename}
        fi
        """

rule assembly_quast:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna"
    output:
        report=f"{OUTPUT_DIR}/cataloging/quast/{{assembly}}/report.tsv"
    params:
        quast_module={QUAST_MODULE},
        outdir=f"{OUTPUT_DIR}/cataloging/quast/{{assembly}}"
    threads: 4
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(4*1024, int(input.size_mb * 5)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(10, int(input.size_mb / 10)) * 2 ** (attempt - 1))
    message: "Calculating assembly statistics for {wildcards.assembly} with QUAST..."
    shell:
        """
        if [ ! -s {input.assembly:q} ]; then
            mkdir -p {params.outdir:q}
            printf 'Assembly\t{wildcards.assembly}\n# contigs\t0\nLargest contig\t0\nTotal length\t0\nGC (%%)\tNA\nN50\t0\nN75\t0\nL50\t0\nL75\t0\n' > {output.report:q}
        else
            module purge
            module load {params.quast_module}
            rm -rf {params.outdir:q}
            quast.py {input.assembly:q} \
                --output-dir {params.outdir:q} \
                --threads {threads} \
                --no-html \
                --no-plots \
                --silent
        fi
        """

rule assembly_map:
    input:
        index=lambda wildcards: f"{OUTPUT_DIR}/cataloging/megahit/{wildcards.assembly}/{wildcards.assembly}.rev.2.bt2",
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        r1=lambda wildcards: PREPROCESSED_TO_READS1[wildcards.sample],
        r2=lambda wildcards: PREPROCESSED_TO_READS2[wildcards.sample]
    output:
        f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}/{{sample}}.bam"
    params:
        bowtie2_module={BOWTIE2_MODULE},
        samtools_module={SAMTOOLS_MODULE},
        basename=lambda wildcards: f"{OUTPUT_DIR}/cataloging/megahit/{wildcards.assembly}/{wildcards.assembly}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 3)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 5)) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.sample} reads to assembly {wildcards.assembly}..."
    shell:
        """
        if [ ! -s {input.assembly} ]; then
            echo "Assembly is empty, skipping mapping..."
            mkdir -p $(dirname {output})
            touch {output}
        else
            module purge
            module load {params.bowtie2_module} {params.samtools_module}
            R1_FILES=$(echo {input.r1} | tr ' ' ',')
            R2_FILES=$(echo {input.r2} | tr ' ' ',')
            bowtie2 -x {params.basename} -1 $R1_FILES -2 $R2_FILES | samtools view -bS - | samtools sort -o {output}
        fi
        """

rule assembly_flagstat:
    input:
        f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}/{{sample}}.bam"
    output:
        f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}/{{sample}}.flagstat.txt"
    params:
        samtools_module={SAMTOOLS_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1*1024, int(input.size_mb / 20)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(20, int(input.size_mb / 100)) * 2 ** (attempt - 1))
    message: "Calculating assembly mapping rate for {wildcards.sample} against {wildcards.assembly}..."
    shell:
        """
        if [ ! -s {input:q} ]; then
            printf '0 + 0 in total (QC-passed reads + QC-failed reads)\n0 + 0 mapped (0.00%% : N/A)\n' > {output:q}
        else
            module purge
            module load {params.samtools_module}
            samtools flagstat -@ {threads} {input:q} > {output:q}
        fi
        """

rule assembly_map_depth:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        bams=lambda wildcards: [
            f"{OUTPUT_DIR}/cataloging/bowtie2/{wildcards.assembly}/{sample}.bam"
            for sample in ASSEMBLY_TO_COVERAGE_SAMPLES[wildcards.assembly]
        ]
    output:
        metabat2=f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}_metabat.depth",
        maxbin2=f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}_maxbin.depth"
    params:
        metabat2_module={METABAT2_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(min(20000,max(15, int(input.size_mb / 100)) * 2 ** (attempt - 1)))
    message: "Calculating mapping states of assembly {wildcards.assembly}..."
    shell:
        """
        if [ ! -s {input.assembly} ]; then
            echo "Assembly is empty, skipping depth calculation..."
            mkdir -p $(dirname {output.metabat2})
            touch {output.metabat2} {output.maxbin2}
        else
            module purge
            module load {params.metabat2_module}
            jgi_summarize_bam_contig_depths --outputDepth {output.metabat2} {input.bams}
            cut -f1,3 {output.metabat2} | tail -n+2 > {output.maxbin2}
        fi
        """

# SemiBin2 and COMEBin keep contig features in memory and read every coverage
# BAM to build the depth profiles used for training, so their footprint is
# driven by the mapping data as much as by the assembly. Sizing on the assembly
# alone under-requests badly for deeply sequenced samples (a 186 MB assembly
# with a 1.8 GB BAM was OOM-killed while training at the previous 16 GB floor).
def _file_size_mb(path):
    try:
        return Path(path).stat().st_size / (1024*1024)
    except OSError:
        return 0


def deep_binner_mem_mb(assembly, bams, attempt):
    """Memory estimate (MB) for the deep-learning binners, scaled by assembly and BAM sizes."""
    bam_mb = [_file_size_mb(bam) for bam in bams] or [0]
    estimate = _file_size_mb(assembly) * 60 + max(bam_mb) * 12 + len(bam_mb) * 1024
    # Scale after the floor so a job that was OOM-killed at the floor actually
    # grows on retry instead of asking for the same memory again.
    return min(1000*1024, max(32*1024, int(estimate)) * 2 ** (attempt - 1))


def deep_binner_runtime(assembly, bams, attempt):
    """Runtime estimate (min) for the deep-learning binners, scaled by assembly and BAM sizes."""
    estimate = _file_size_mb(assembly) * 2 + sum(_file_size_mb(bam) for bam in bams) / 10
    return min(20000, max(120, int(estimate)) * 2 ** (attempt - 1))


def binning_skip_reason(assembly):
    """Why an assembly should not be binned, or an empty string when it should.

    The binners share one policy, evaluated per job rather than at DAG time
    because the assembly does not exist yet when the DAG is built.
    """
    size_mb = _file_size_mb(assembly)
    if size_mb <= 0:
        return "Assembly is empty"
    if size_mb < MIN_BINNING_ASSEMBLY_MB:
        return f"Assembly is smaller than {MIN_BINNING_ASSEMBLY_MB} MB"
    return ""


rule metabat2:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        depth=f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}_metabat.depth"
    output:
        f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.tsv"
    params:
        metabat2_module={METABAT2_MODULE},
        raw_cls=f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.raw.tsv",
        skip_reason=lambda wildcards, input: binning_skip_reason(input.assembly)
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 50)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 5)) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using metabat2..."
    shell:
        """
        if [ -n "{params.skip_reason}" ]; then
            echo "{params.skip_reason}, skipping metabat2..."
            mkdir -p $(dirname {output})
            touch {output}
        else
            module purge
            module load {params.metabat2_module}
            # --saveCls writes one row per contig, including every contig
            # metabat2 left unbinned (cluster 0) and every contig below -m 1500.
            # Handing those to binette would turn them into a single junk bin
            # that drags the whole assembly through gene prediction and
            # diamond, so keep only the contigs that were actually binned.
            # The header filter covers metabat2 > 2.17, which prepends a
            # "ContigName ClusterId" line to the membership matrix.
            metabat2 -i {input.assembly} -a {input.depth} -o {params.raw_cls} -m 1500 --saveCls --noBinOut
            awk -F'\t' '$1 != "ContigName" && $2 != "0"' {params.raw_cls} > {output}
            rm -f {params.raw_cls}
        fi
        """

rule maxbin2:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        depth=f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}_maxbin.depth"
    output:
        f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.summary"
    params:
        maxbin2_module={MAXBIN2_MODULE},
        hmmer_module={HMMER_MODULE},
        basename=f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}",
        skip_reason=lambda wildcards, input: binning_skip_reason(input.assembly)
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 50)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 3)) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using maxbin2..."
    shell:
        """
        if [ -n "{params.skip_reason}" ]; then
            echo "{params.skip_reason}, skipping maxbin2..."
            mkdir -p $(dirname {output})
            touch {output}
        else
            MODULEPATH=/opt/shared_software/shared_envmodules/modules:$MODULEPATH \
            module purge
            module load {params.maxbin2_module} {params.hmmer_module}
            rm -rf {params.basename}*
            run_MaxBin.pl -contig {input.assembly} -abund {input.depth} -max_iteration 10 -out {params.basename} -min_contig_length 1500
        fi
        """

rule maxbin2_table:
    input:
        f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.summary"
    output:
        f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}/{{assembly}}.tsv"
    params:
        package_dir={PACKAGE_DIR},
        fastadir=f"{OUTPUT_DIR}/cataloging/maxbin2/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 10)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 5)) * 2 ** (attempt - 1))
    shell:
        """
        if [ ! -s {input} ]; then
            touch {output}
        else
            python {params.package_dir}/workflow/scripts/fastas_to_bintable.py -d {params.fastadir} -e fasta -o {output}
        fi
        """


rule semibin2:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        bam=lambda wildcards: [
            f"{OUTPUT_DIR}/cataloging/bowtie2/{wildcards.assembly}/{sample}.bam"
            for sample in ASSEMBLY_TO_COVERAGE_SAMPLES[wildcards.assembly]
            ]
    output:
        f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/contig_bins.tsv"
    params:
        outdir=f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}",
        skip_reason=lambda wildcards, input: binning_skip_reason(input.assembly)
    threads: 8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/semibin.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(deep_binner_mem_mb(input.assembly, input.bam, attempt)),
        runtime=lambda wildcards, input, attempt: cap_runtime(deep_binner_runtime(input.assembly, input.bam, attempt)),
        slurm_partition="gpuqueue",
        slurm_extra="--gres=gpu:1"
    message: "Binning contigs from assembly {wildcards.assembly} using semibin2..."
    shell:
        """
        # Prevent the snakemake module's PYTHONPATH (and ~/.local) from shadowing
        # the conda env's own Python packages, e.g. an older narwhals that
        # scikit-learn >= 1.9 cannot import.
        unset PYTHONPATH
        export PYTHONNOUSERSITE=1
        if [ -n "{params.skip_reason}" ]; then
            echo "{params.skip_reason}, skipping semibin2..."
            mkdir -p {params.outdir}
            touch {output}
        else
            SemiBin2 single_easy_bin -i {input.assembly} -b {input.bam} -o {params.outdir} -m 1500 -t {threads} --compression none
        fi
        """

rule semibin2_table:
    input:
        f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/contig_bins.tsv"
    output:
        f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/{{assembly}}.tsv"
    params:
        package_dir={PACKAGE_DIR},
        fastadir=f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}/output_bins"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1*1024, int(input.size_mb * 10)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(5, int(input.size_mb / 5)) * 2 ** (attempt - 1))
    shell:
        """
        if [ ! -s {input} ]; then
            touch {output}
        else
            python {params.package_dir}/workflow/scripts/fastas_to_bintable.py -d {params.fastadir} -e fa -o {output}
        fi
        """

rule comebin:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        bam=lambda wildcards: [
            f"{OUTPUT_DIR}/cataloging/bowtie2/{wildcards.assembly}/{sample}.bam"
            for sample in ASSEMBLY_TO_COVERAGE_SAMPLES[wildcards.assembly]
            ]
    output:
        f"{OUTPUT_DIR}/cataloging/comebin/{{assembly}}/comebin_res/comebin_res.tsv"
    params:
        bamdir=f"{OUTPUT_DIR}/cataloging/comebin/{{assembly}}_bams",
        outdir=f"{OUTPUT_DIR}/cataloging/comebin/{{assembly}}",
        skip_reason=lambda wildcards, input: binning_skip_reason(input.assembly)
    threads: 8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/comebin.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(deep_binner_mem_mb(input.assembly, input.bam, attempt)),
        runtime=lambda wildcards, input, attempt: cap_runtime(deep_binner_runtime(input.assembly, input.bam, attempt)),
        slurm_partition="gpuqueue",
        slurm_extra="--gres=gpu:1"
    message: "Binning contigs from assembly {wildcards.assembly} using comebin..."
    shell:
        """
        # Prevent the snakemake module's PYTHONPATH (and ~/.local) from shadowing
        # the conda env's own Python packages.
        unset PYTHONPATH
        unset PYTHONHOME
        export PYTHONNOUSERSITE=1
        # run_comebin.sh calls a bare `python`. When a conda env is already
        # active in the submitting shell, `conda activate` swaps the new env's
        # bin into the *position* of the old one instead of prepending it, so
        # the snakemake module's Python (3.11, no torch) stays ahead of this
        # env's Python 3.7 in PATH. Put the env's bin first explicitly.
        if [ -n "${{CONDA_PREFIX:-}}" ]; then
            export PATH="$CONDA_PREFIX/bin:$PATH"
        fi
        # COMEBin seeds its clustering with single-copy marker genes found by
        # FragGeneScan and HMMsearch. Assemblies too small to carry those markers
        # make that stage abort ("produced no proteins", "zero marker hits",
        # "empty seed set"), which used to kill the whole workflow. Skip them up
        # front instead of spending GPU hours on a run that cannot finish.
        if [ -n "{params.skip_reason}" ]; then
            echo "{params.skip_reason}, skipping comebin..."
            mkdir -p $(dirname {output})
            touch {output}
        else
            rm -rf {params.outdir} {params.bamdir}
            mkdir -p {params.bamdir} {params.outdir}
            for bam in {input.bam}; do
                ln -s "$(realpath "$bam")" "{params.bamdir}/$(basename "$bam")"
            done

            # COMEBin writes its marker byproducts next to the assembly and reuses
            # them without checking that they hold anything. An empty one left by a
            # previous failed attempt would make every retry fail the same way, so
            # drop the empties and let COMEBin regenerate them.
            COMEBIN_SEED="{input.assembly}.bacar_marker.2quarter_lencutoff_1001.seed"
            for byproduct in "{input.assembly}.frag.faa" "{input.assembly}.bacar_marker.hmmout" "$COMEBIN_SEED"; do
                if [ -f "$byproduct" ] && [ ! -s "$byproduct" ]; then
                    rm -f "$byproduct"
                fi
            done

            COMEBIN_LOG="{params.outdir}/comebin.log"
            set +e
            run_comebin.sh \
                -a {input.assembly} \
                -p {params.bamdir} \
                -t {threads} \
                -o {params.outdir} 2>&1 | tee "$COMEBIN_LOG"
            COMEBIN_STATUS=${{PIPESTATUS[0]}}
            set -e

            # The size threshold cannot catch every marker-poor assembly: a large
            # but heavily fragmented one can still carry no single-copy markers.
            # COMEBin signals that either by aborting outright (1.1.0) or by
            # logging the failed marker step and walking away without a result
            # (1.0.4, which exits 0 and leaves snakemake to die on the missing
            # output). Detect both, and treat only that specific failure as
            # "nothing to bin": every other failure (OOM, GPU, walltime) still
            # propagates so snakemake can retry it with more resources.
            MARKER_FAILURE=0
            if grep -qE "FragGeneScan failed|Hmmsearch failed|markerCmd failed|produced no proteins|zero marker hits|empty seed set|no usable seeds" "$COMEBIN_LOG"; then
                MARKER_FAILURE=1
            fi
            if [ -f "$COMEBIN_SEED" ] && [ ! -s "$COMEBIN_SEED" ]; then
                MARKER_FAILURE=1
            fi

            if [ "$COMEBIN_STATUS" -ne 0 ] || [ ! -f {output} ]; then
                if [ "$MARKER_FAILURE" -eq 1 ]; then
                    echo "COMEBin found no marker genes in assembly {wildcards.assembly}, exporting an empty comebin_res.tsv..."
                    mkdir -p $(dirname {output})
                    touch {output}
                    exit 0
                fi
                if [ "$COMEBIN_STATUS" -ne 0 ]; then
                    exit "$COMEBIN_STATUS"
                fi
                echo "COMEBin reported success but wrote no {output}. See $COMEBIN_LOG." >&2
                exit 1
            fi
        fi
        """

rule comebin_table:
    input:
        f"{OUTPUT_DIR}/cataloging/comebin/{{assembly}}/comebin_res/comebin_res.tsv"
    output:
        f"{OUTPUT_DIR}/cataloging/comebin/{{assembly}}/{{assembly}}.tsv"
    params:
        package_dir={PACKAGE_DIR},
        fastadir=f"{OUTPUT_DIR}/cataloging/comebin/{{assembly}}/comebin_res/comebin_res_bins"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1*1024, int(input.size_mb * 10)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(5, int(input.size_mb / 5)) * 2 ** (attempt - 1))
    shell:
        """
        if [ ! -s {input} ]; then
            touch {output}
        else
            python {params.package_dir}/workflow/scripts/fastas_to_bintable.py -d {params.fastadir} -e fa -o {output}
        fi
        """

# Script to calculate resources based on the number of bins
_row_count_cache = {}
def row_count(path):
    """Return number of data rows (excluding header) in a TSV, caching the result."""
    if path not in _row_count_cache:
        with open(path) as f:
            _row_count_cache[path] = max(0, sum(1 for _ in f))
    return _row_count_cache[path]


def row_count_sum(paths):
    return sum(row_count(path) for path in paths)


checkpoint binette:
    input:
        binner_tables=selected_binner_tables,
        fasta=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna"
    output:
        f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins_quality_reports.tsv"
    params:
        checkm_db = {CHECKM2_DB},
        diamond_module = {DIAMOND_MODULE},
        checkm2_module = {CHECKM2_MODULE},
        binette_module = {BINETTE_MODULE},
        min_completeness={MIN_COMPLETENESS},
        max_contamination={MAX_CONTAMINATION},
        min_bin_length={MIN_BIN_LENGTH},
        max_bin_length={MAX_BIN_LENGTH},
        bin_prefix={BINETTE_BIN_PREFIX},
        report_header={BINETTE_REPORT_HEADER},
        outdir=f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}"
    threads: 8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/cataloging.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(min(1000*1024,max(32*1024, row_count_sum(input.binner_tables)) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(min(20000,max(15, int(input.size_mb)) * 2 ** (attempt - 1)))
    message: "Refining bins from assembly {wildcards.assembly} using binette..."
    shell:
        """
        # Remove empty input files from the list
        VALID_TSV_FILES=()
        for table in {input.binner_tables:q}; do
            if [ -s "$table" ]; then
                VALID_TSV_FILES+=("$table")
            fi
        done

        mkdir -p {params.outdir}

        # Ensure at least one valid TSV file exists
        if [ "${{#VALID_TSV_FILES[@]}}" -eq 0 ]; then
            echo "No valid TSV input files for binette, skipping..."
            printf "{params.report_header}" > {output}
            exit 0
        fi

        # Run binette only with non-empty TSV files
        DIAMOND_RESULT_TSV="{params.outdir}/temporary_files/diamond_result.tsv"
        DIAMOND_RESULT_TSV_GZ="{params.outdir}/temporary_files/diamond_result.tsv.gz"
        EMPTY_OUTPUT_HEADER='{params.report_header}'

        set +e
        binette --contig2bin_tables "${{VALID_TSV_FILES[@]}}" \
                --contigs {input.fasta} \
                --outdir {params.outdir} \
                --checkm2_db {params.checkm_db} \
                --prefix {params.bin_prefix} \
                --min_completeness {params.min_completeness} \
                --max_contamination {params.max_contamination} \
                --min_length {params.min_bin_length} \
                --max_length {params.max_bin_length} \
                --threads {threads}
        BINETTE_STATUS=$?
        set -e

        if [ "$BINETTE_STATUS" -eq 0 ]; then
            exit 0
        fi

        EMPTY_DIAMOND_RESULT=0
        if [ -f "$DIAMOND_RESULT_TSV" ] && [ ! -s "$DIAMOND_RESULT_TSV" ]; then
            EMPTY_DIAMOND_RESULT=1
        elif [ -f "$DIAMOND_RESULT_TSV_GZ" ] && [ "$(gzip -cd "$DIAMOND_RESULT_TSV_GZ" 2>/dev/null | wc -c | tr -d '[:space:]')" = "0" ]; then
            EMPTY_DIAMOND_RESULT=1
        fi

        if [ "$EMPTY_DIAMOND_RESULT" -eq 1 ]; then
            echo "Binette produced an empty diamond_result file, exporting an empty final_bins_quality_reports.tsv..."
            printf "$EMPTY_OUTPUT_HEADER" > {output}
            exit 0
        fi

        exit "$BINETTE_STATUS"
        """

# Regenerate the bin_id wildcard based on the checkpoint results
def get_bin_fna_sep(wildcards):
    checkpoint_output = checkpoints.binette.get(**wildcards).output[0]
    cluster_ids = get_bin_ids_from_tsv(checkpoint_output)
    return f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins/{BINETTE_BIN_PREFIX}_bin{wildcards.bin_id}.fa"

rule rename_bins:
    input:
        lambda wildcards: get_bin_fna_sep(wildcards)
    output:
        f"{OUTPUT_DIR}/cataloging/final/{{assembly}}/{{assembly}}_bin_{{bin_id}}.fa"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1*1024, int(input.size_mb * 10)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(2, int(input.size_mb)) * 2 ** (attempt - 1))
    message: "Copying bin {wildcards.bin_id} from assembly {wildcards.assembly}..."
    shell:
        """
        cp {input} {output}
        """

rule move_metadata:
    input:
        f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins_quality_reports.tsv"
    output:
        f"{OUTPUT_DIR}/cataloging/final/{{assembly}}.tsv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(8*1024 * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(10 * 2 ** (attempt - 1))
    message: "Exporting bin metadata from assembly {wildcards.assembly}..."
    shell:
        """
        cp {input} {output}
        """

# The binette bins keep the contig names assigned by the `assembly` rule
# (`<assembly>_<n>`), so the table below is the explicit link between every
# final bin and the assembly contigs it is made of.
def get_final_bin_fastas(wildcards):
    checkpoint_output = checkpoints.binette.get(**wildcards).output[0]
    bin_ids = get_bin_ids_from_tsv(checkpoint_output)
    return expand(
        f"{OUTPUT_DIR}/cataloging/final/{{assembly}}/{{assembly}}_bin_{{bin_id}}.fa",
        assembly=wildcards.assembly,
        bin_id=bin_ids
    )

rule contig_to_bin:
    input:
        bins=get_final_bin_fastas
    output:
        f"{OUTPUT_DIR}/cataloging/contig_to_bin/{{assembly}}.tsv"
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1*1024, int(input.size_mb)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(5, int(input.size_mb / 10)) * 2 ** (attempt - 1))
    message: "Mapping contigs to bins from assembly {wildcards.assembly}..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/contig_to_bin.py \
            --assembly {wildcards.assembly:q} \
            -o {output:q} \
            {input.bins:q}
        """

rule all_bins:
    input:
        expand(f"{OUTPUT_DIR}/cataloging/final/{{assembly}}.tsv", assembly=assemblies)
    output:
        paths=f"{OUTPUT_DIR}/cataloging/final/all_bin_paths.txt",
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv"
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1*1024, int(input.size_mb * 10)) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(2, int(input.size_mb)) * 2 ** (attempt - 1))
    message: "Generating bin path file..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/all_bin_paths.py {input} -o {output.paths}
        python {params.package_dir}/workflow/scripts/all_bin_metadata.py {input} -o {output.metadata}
        """

rule all_contig_to_bin:
    input:
        expand(f"{OUTPUT_DIR}/cataloging/contig_to_bin/{{assembly}}.tsv", assembly=assemblies)
    output:
        f"{OUTPUT_DIR}/cataloging/final/all_contig_to_bin.csv"
    localrule: True
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(1*1024 * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(10 * 2 ** (attempt - 1))
    message: "Combining contig-to-bin tables..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/contig_to_bin.py \
            --merge \
            -o {output:q} \
            {input:q}
        """

rule cataloging_stats:
    input:
        assembly_to_samples=f"{OUTPUT_DIR}/data/assembly_to_samples.json",
        quast=expand(f"{OUTPUT_DIR}/cataloging/quast/{{assembly}}/report.tsv", assembly=assemblies),
        flagstats=[
            f"{OUTPUT_DIR}/cataloging/bowtie2/{assembly}/{sample}.flagstat.txt"
            for assembly, samples in ASSEMBLY_TO_COVERAGE_SAMPLES.items()
            for sample in samples
        ],
        metabat2=selected_binner_expand("metabat"),
        maxbin2=selected_binner_expand("maxbin"),
        semibin2=selected_binner_expand("semibin"),
        comebin=selected_binner_expand("comebin"),
        bins=expand(f"{OUTPUT_DIR}/cataloging/final/{{assembly}}.tsv", assembly=assemblies)
    output:
        f"{OUTPUT_DIR}/cataloging.tsv"
    localrule: True
    params:
        package_dir={PACKAGE_DIR},
        binette_report_root=f"{OUTPUT_DIR}/cataloging/binette"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(1*1024 * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: cap_runtime(5 * 2 ** (attempt - 1))
    message: "Creating cataloging stats..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/cataloging_stats.py \
            --assembly-to-samples {input.assembly_to_samples:q} \
            --quast {input.quast} \
            --flagstat {input.flagstats} \
            --metabat2 {input.metabat2} \
            --maxbin2 {input.maxbin2} \
            --semibin2 {input.semibin2} \
            --comebin {input.comebin} \
            --binette-report-root {params.binette_report_root:q} \
            --bins {input.bins} \
            -o {output:q}
        """
