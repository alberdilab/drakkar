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
BEDTOOLS_MODULE = config["BEDTOOLS_MODULE"]
HMMER_MODULE = config["HMMER_MODULE"]
COMEBIN_MODULE = config["COMEBIN_MODULE"]
SEMIBIN2_MODULE = config["SEMIBIN2_MODULE"]
DIAMOND_MODULE = config["DIAMOND_MODULE"]
CHECKM2_MODULE = config["CHECKM2_MODULE"]
BINETTE_MODULE = config["BINETTE_MODULE"]

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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(min(1020*1024,max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 20) * 2 ** (attempt - 1)))
    message: "Assembling {wildcards.assembly}..."
    shell:
        """
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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 5) * 2 ** (attempt - 1)))
    message: "Indexing assembly {wildcards.assembly}..."
    shell:
        """
        if [ ! -s {input} ]; then
            echo "Assembly is empty, skipping bowtie2-build..."
            touch {output.index}
        else
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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(4*1024, int(input.size_mb * 5) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(10, int(input.size_mb / 10) * 2 ** (attempt - 1)))
    message: "Calculating assembly statistics for {wildcards.assembly} with QUAST..."
    shell:
        """
        if [ ! -s {input.assembly:q} ]; then
            mkdir -p {params.outdir:q}
            printf 'Assembly\t{wildcards.assembly}\n# contigs\t0\nLargest contig\t0\nTotal length\t0\nGC (%%)\tNA\nN50\t0\nN75\t0\nL50\t0\nL75\t0\n' > {output.report:q}
        else
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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 3) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 5) * 2 ** (attempt - 1)))
    message: "Mapping {wildcards.sample} reads to assembly {wildcards.assembly}..."
    shell:
        """
        if [ ! -s {input.assembly} ]; then
            echo "Assembly is empty, skipping mapping..."
            mkdir -p $(dirname {output})
            touch {output}
        else
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
        mem_mb=cap_mem_mb(1*1024),
        runtime=cap_runtime(5)
    message: "Calculating assembly mapping rate for {wildcards.sample} against {wildcards.assembly}..."
    shell:
        """
        if [ ! -s {input:q} ]; then
            printf '0 + 0 in total (QC-passed reads + QC-failed reads)\n0 + 0 mapped (0.00%% : N/A)\n' > {output:q}
        else
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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(min(20000,max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))))
    message: "Calculating mapping states of assembly {wildcards.assembly}..."
    shell:
        """
        if [ ! -s {input.assembly} ]; then
            echo "Assembly is empty, skipping depth calculation..."
            mkdir -p $(dirname {output.metabat2})
            touch {output.metabat2} {output.maxbin2}
        else
            module load {params.metabat2_module}
            jgi_summarize_bam_contig_depths --outputDepth {output.metabat2} {input.bams}
            cut -f1,3 {output.metabat2} | tail -n+2 > {output.maxbin2}
        fi
        """

rule metabat2:
    input:
        assembly=f"{OUTPUT_DIR}/cataloging/megahit/{{assembly}}/{{assembly}}.fna",
        depth=f"{OUTPUT_DIR}/cataloging/bowtie2/{{assembly}}_metabat.depth"
    output:
        f"{OUTPUT_DIR}/cataloging/metabat2/{{assembly}}/{{assembly}}.tsv"
    params:
        metabat2_module={METABAT2_MODULE}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 5) * 2 ** (attempt - 1)))
    message: "Binning contigs from assembly {wildcards.assembly} using metabat2..."
    shell:
        """
        if [ ! -s {input.assembly} ]; then
            echo "Assembly is empty, skipping metabat2..."
            mkdir -p $(dirname {output})
            touch {output}
        else
            module load {params.metabat2_module}
            metabat2 -i {input.assembly} -a {input.depth} -o {output} -m 1500 --saveCls --noBinOut
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
        assembly_size_mb=lambda wildcards, input: int(Path(input.assembly).stat().st_size / (1024*1024))
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 3) * 2 ** (attempt - 1)))
    message: "Binning contigs from assembly {wildcards.assembly} using maxbin2..."
    shell:
        """
        if (( {params.assembly_size_mb} < 10 )); then
            echo "Assembly is smaller than 10 MB, skipping maxbin2..."
            mkdir -p $(dirname {output})
            touch {output}
        else
            MODULEPATH=/opt/shared_software/shared_envmodules/modules:$MODULEPATH \
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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(15, int(input.size_mb / 5) * 2 ** (attempt - 1)))
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
        semibin2_module={SEMIBIN2_MODULE},
        hmmer_module={HMMER_MODULE},
        bedtools_module={BEDTOOLS_MODULE},
        outdir=f"{OUTPUT_DIR}/cataloging/semibin2/{{assembly}}",
        assembly_size_mb=lambda wildcards, input: int(Path(input.assembly).stat().st_size / (1024*1024))
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(min(1000*1024,max(8*1024, int(input.size_mb * 30) * 2 ** (attempt - 1)))),
        runtime=lambda wildcards, input, attempt: cap_runtime(min(20000,max(15, int(input.size_mb / 2) * 2 ** (attempt - 1))))
    message: "Binning contigs from assembly {wildcards.assembly} using semibin2..."
    shell:
        """
        if (( {params.assembly_size_mb} < 10 )); then
            echo "Assembly is smaller than 10 MB, skipping semibin2..."
            mkdir -p {params.outdir}
            touch {output}
        else
            module load {params.semibin2_module} {params.bedtools_module} {params.hmmer_module}
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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(5, int(input.size_mb / 5) * 2 ** (attempt - 1)))
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
        comebin_module={COMEBIN_MODULE},
        bamdir=f"{OUTPUT_DIR}/cataloging/comebin/{{assembly}}_bams",
        outdir=f"{OUTPUT_DIR}/cataloging/comebin/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(min(1000*1024,max(8*1024, int(input.size_mb * 30) * 2 ** (attempt - 1)))),
        runtime=lambda wildcards, input, attempt: cap_runtime(min(20000,max(15, int(input.size_mb / 2) * 2 ** (attempt - 1))))
    message: "Binning contigs from assembly {wildcards.assembly} using comebin..."
    shell:
        """
        if [ ! -s {input.assembly} ]; then
            echo "Assembly is empty, skipping comebin..."
            mkdir -p $(dirname {output})
            touch {output}
        else
            module load {params.comebin_module}
            rm -rf {params.outdir} {params.bamdir}
            mkdir -p {params.bamdir}
            for bam in {input.bam}; do
                ln -s "$(realpath "$bam")" "{params.bamdir}/$(basename "$bam")"
            done
            run_comebin.sh \
                -a {input.assembly} \
                -p {params.bamdir} \
                -t {threads} \
                -o {params.outdir}
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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(5, int(input.size_mb / 5) * 2 ** (attempt - 1)))
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
        outdir=f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}"
    threads: 8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/cataloging.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(min(1000*1024,max(32*1024, row_count_sum(input.binner_tables) * 2 ** (attempt - 1)))),
        runtime=lambda wildcards, input, attempt: cap_runtime(min(20000,max(15, int(input.size_mb) * 2 ** (attempt - 1))))
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
            printf "bin_id\tcompleteness\tcontamination\tscore\tsize\tN50\tcontig_count\n" > {output}
            exit 0
        fi

        # Run binette only with non-empty TSV files
        DIAMOND_RESULT_TSV="{params.outdir}/temporary_files/diamond_result.tsv"
        DIAMOND_RESULT_TSV_GZ="{params.outdir}/temporary_files/diamond_result.tsv.gz"
        EMPTY_OUTPUT_HEADER='bin_id\tcompleteness\tcontamination\tscore\tsize\tN50\tcontig_count\n'

        set +e
        binette --contig2bin_tables "${{VALID_TSV_FILES[@]}}" \
                --contigs {input.fasta} \
                --outdir {params.outdir} \
                --checkm2_db {params.checkm_db} \
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
    return f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins/bin_{wildcards.bin_id}.fa"

rule rename_bins:
    input:
        lambda wildcards: get_bin_fna_sep(wildcards)
    output:
        f"{OUTPUT_DIR}/cataloging/final/{{assembly}}/{{assembly}}_bin_{{bin_id}}.fa"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(2, int(input.size_mb) * 2 ** (attempt - 1)))
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
        mem_mb=cap_mem_mb(8*1024),
        runtime=cap_runtime(10)
    message: "Exporting bin metadata from assembly {wildcards.assembly}..."
    shell:
        """
        cp {input} {output}
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
        mem_mb=lambda wildcards, input, attempt: cap_mem_mb(max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: cap_runtime(max(2, int(input.size_mb) * 2 ** (attempt - 1)))
    message: "Generating bin path file..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/all_bin_paths.py {input} -o {output.paths}
        python {params.package_dir}/workflow/scripts/all_bin_metadata.py {input} -o {output.metadata}
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
        mem_mb=cap_mem_mb(1*1024),
        runtime=cap_runtime(5)
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
