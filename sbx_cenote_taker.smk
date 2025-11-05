try:
    SBX_CENOTE_TAKER_VERSION = get_ext_version("sbx_cenote_taker")
except (NameError, ValueError):
    # For backwards compatibility with older versions of Sunbeam
    SBX_CENOTE_TAKER_VERSION = "0.0.0"
VIRUS_FP = output_subdir(Cfg, "virus")


def get_extension_path() -> Path:
    return Path(__file__).parent.resolve()


def cenote_output() -> Path:
    return VIRUS_FP / "cenote_taker" / "{sample}.fasta"


rule all_cenote_taker:
    input:
        expand(
            VIRUS_FP / "alignments" / "{sample}.gene_coverage.tsv",
            sample=Samples.keys(),
        ),
        expand(
            VIRUS_FP / "blastx" / "{sample}.btf",
            sample=Samples.keys(),
        ),
        VIRUS_FP / "summary" / "all_align_summary.txt",


rule cenote_taker:
    input:
        contigs=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
    output:
        contigs=VIRUS_FP / "cenote_taker" / "{sample}" / "final.contigs.fasta",
        summary=VIRUS_FP
        / "cenote_taker"
        / "{sample}"
        / "{sample}"
        / "{sample}_CONTIG_SUMMARY.tsv",
    benchmark:
        BENCHMARK_FP / "cenote_taker_{sample}.tsv"
    log:
        LOG_FP / "cenote_taker_{sample}.log",
    params:
        out_dir=str(VIRUS_FP / "cenote_taker"),
        sample="{sample}",
        db_fp=Cfg["sbx_cenote_taker"]["cenote_taker_db"],
    resources:
        mem_mb=24000,
        runtime=720,
    conda:
        "envs/cenote_taker_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-cenote-taker"
    shell:
        """
        SAMPLE={params.sample}
        if [[ ${{#SAMPLE}} -lt 18 ]] && [[ {params.sample} =~ ^[a-zA-Z0-9_]+$ ]]; then
            echo "Sample name format is valid" >> {log}
        else
            echo "Cenote-Taker requires a sample name that is less than 18 characters and contains only alphanumeric characters and underscores" >> {log}
            exit 1
        fi

        if [ -s {input.contigs} ]; then
            echo "Contigs file exists and is not empty" >> {log}
        else
            echo "Contigs file is empty" >> {log}
            touch {output.contigs} {output.summary}
            exit 0
        fi

        if [ ! -d {params.db_fp} ] || [ ! "$(ls -A {params.db_fp})" ]; then
            echo "Cenote-Taker database path {params.db_fp} is missing or empty" >> {log}
            exit 1
        fi

        cd {params.out_dir}
        export CENOTE_DBS={params.db_fp}
        cenotetaker3 --contigs {input.contigs} -r {params.sample} -p T --lin_minimum_hallmark_genes 2 >> {log} 2>&1
        """


rule filter_cenote_contigs:
    input:
        contigs=VIRUS_FP / "cenote_taker" / "{sample}" / "final.contigs.fasta",
        summary=VIRUS_FP
        / "cenote_taker"
        / "{sample}"
        / "{sample}"
        / "{sample}_CONTIG_SUMMARY.tsv",
    output:
        VIRUS_FP / "cenote_taker" / "{sample}.fasta",
    params:
        include_phages=Cfg["sbx_cenote_taker"]["include_phages"],
    script:
        "scripts/filter_cenote_contigs.py"


rule build_virus_index:
    input:
        cenote_output(),
    output:
        str(cenote_output()) + ".1.bt2",  # Don't use f-string, broken with python 3.12
    conda:
        "envs/sbx_cenote_taker.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-sbx-cenote-taker"
    threads: Cfg["sbx_cenote_taker"]["bowtie2_build_threads"]
    shell:
        "bowtie2-build --threads {threads} -f {input} {input}"


rule align_virus_reads:
    input:
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
        index=str(cenote_output()) + ".1.bt2",  # Don't use f-string, broken with python 3.12
    output:
        temp(VIRUS_FP / "alignments" / "{sample}.sam"),
    params:
        index=str(cenote_output()),
    threads: 6
    conda:
        "envs/sbx_cenote_taker.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-sbx-cenote-taker"
    shell:
        "bowtie2 -q --local -t --very-sensitive-local --threads {threads} --no-mixed --no-discordant -x {params.index} -1 {input.r1} -2 {input.r2} -S {output}"


rule process_virus_alignment:
    input:
        VIRUS_FP / "alignments" / "{sample}.sam",
    output:
        bam=temp(VIRUS_FP / "alignments" / "{sample}.bam"),
        sorted=temp(VIRUS_FP / "alignments" / "{sample}.sorted.bam"),
        bai=temp(VIRUS_FP / "alignments" / "{sample}.sorted.bam.bai"),
    params:
        target=str(cenote_output()),
    conda:
        "envs/sbx_cenote_taker.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-sbx-cenote-taker"
    shell:
        """
        samtools view -bT {params.target} {input} > {output.bam}
        samtools sort -o {output.sorted} {output.bam}
        samtools index {output.sorted} {output.bai}
        """


rule calculate_mapping_stats:
    input:
        bam=VIRUS_FP / "alignments" / "{sample}.sorted.bam",
        idx=VIRUS_FP / "alignments" / "{sample}.sorted.bam.bai",
    output:
        VIRUS_FP / "alignments" / "{sample}.sorted.idxstats.tsv",
    conda:
        "envs/sbx_cenote_taker.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-sbx-cenote-taker"
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """


rule virus_mpileup:
    input:
        bam=VIRUS_FP / "alignments" / "{sample}.sorted.bam",
        idx=VIRUS_FP / "alignments" / "{sample}.sorted.bam.bai",
        contigs=cenote_output(),
    output:
        VIRUS_FP / "alignments" / "{sample}.mpileup",
    conda:
        "envs/sbx_cenote_taker.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-sbx-cenote-taker"
    shell:
        """
        samtools mpileup -f {input.contigs} {input.bam} > {output}
        """


rule filter_virus_coverage:
    input:
        fa=cenote_output(),
        idx=VIRUS_FP / "alignments" / "{sample}.sorted.idxstats.tsv",
    output:
        VIRUS_FP / "final_{sample}_contigs.fasta",
    log:
        LOG_FP / "filter_virus_coverage_{sample}.log",
    script:
        "scripts/filter_virus_coverage.py"


rule virus_blastx:
    """Run blastx on untranslated genes against a target db and write to blast tabular format."""
    input:
        VIRUS_FP / "final_{sample}_contigs.fasta",
    output:
        VIRUS_FP / "blastx" / "{sample}.btf",
    benchmark:
        BENCHMARK_FP / "run_virus_blastx_{sample}.tsv"
    log:
        LOG_FP / "run_virus_blastx_{sample}.log",
    params:
        blast_db=Cfg["sbx_cenote_taker"]["blast_db"],
    threads: Cfg["sbx_cenote_taker"]["blastx_threads"]
    resources:
        mem_mb=24000,
        runtime=720,
    conda:
        "envs/sbx_cenote_taker.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-sbx-cenote-taker"
    shell:
        """
        if [ -s {input} ]; then
            export BLASTDB=$(dirname {params.blast_db})
            blastx \
            -query {input} \
            -db $(basename {params.blast_db}) \
            -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
            -num_threads {threads} \
            -evalue 0.05 \
            -max_target_seqs 100 \
            -out {output} \
            2>&1 | tee {log}
        else
            echo "Caught empty query" >> {log}
            touch {output}
        fi
        """


rule calculate_coverage:
    input:
        bam=VIRUS_FP / "alignments" / "{sample}.sorted.bam",
        idx=VIRUS_FP / "alignments" / "{sample}.sorted.bam.bai",
    output:
        VIRUS_FP / "alignments" / "{sample}.genomecoverage.txt",
    params:
        ext_fp=str(get_extension_path()),
    conda:
        "envs/sbx_cenote_taker.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-sbx-cenote-taker"
    shell:
        """
        samtools view -b {input.bam} | genomeCoverageBed -ibam stdin | grep -v 'genome'| perl {params.ext_fp}/scripts/coverage_counter.pl > {output}
        """


rule combine_coverage_stats:
    input:
        cov=VIRUS_FP / "alignments" / "{sample}.genomecoverage.txt",
        stats=VIRUS_FP / "alignments" / "{sample}.sorted.idxstats.tsv",
    output:
        VIRUS_FP / "alignments" / "{sample}.align.summary.txt",
    benchmark:
        BENCHMARK_FP / "combine_coverage_stats_{sample}.tsv"
    log:
        LOG_FP / "combine_coverage_stats_{sample}.log",
    params:
        ext_fp=str(get_extension_path()),
    conda:
        "envs/r_env.yml"
    container:
        "docker://r-base:latest"
    shell:
        """
        Rscript {params.ext_fp}/scripts/combine_coverage_stats.R {input.cov} {input.stats} {output} 2>&1 | tee {log}
        """


rule virus_coverage_per_gene:
    input:
        mpileup=VIRUS_FP / "alignments" / "{sample}.mpileup",
        btf=VIRUS_FP / "blastx" / "{sample}.btf",
    output:
        tsv=VIRUS_FP / "alignments" / "{sample}.gene_coverage.tsv",
    params:
        contigs=cenote_output(),
    conda:
        "envs/sbx_cenote_taker.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-sbx-cenote-taker"
    script:
        "scripts/virus_coverage_per_gene.py"


rule all_summary:
    input:
        expand(
            VIRUS_FP / "alignments" / "{sample}.align.summary.txt",
            sample=Samples.keys(),
        ),
    output:
        VIRUS_FP / "summary" / "all_align_summary.txt",
    shell:
        """
        echo -e "Sample\tAlignTarget\tFractionCoverage\tTargetLength\tMappedReads" > {output}
        cat {input} >> {output}
        """
