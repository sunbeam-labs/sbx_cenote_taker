try:
    SBX_CENOTE_TAKER_VERSION = get_ext_version("sbx_cenote_taker")
except (NameError, ValueError):
    # For backwards compatibility with older versions of Sunbeam
    SBX_CENOTE_TAKER_VERSION = "0.0.0"
VIRUS_FP = output_subdir(Cfg, "virus")


def get_extension_path() -> Path:
    return Path(__file__).parent.parent.resolve()
def get_rules_path() -> Path:
    return Path(__file__).resolve()


rule all_cenote_taker:
    input:
        expand(
            VIRUS_FP / "cenote_taker" / "{sample}" / "{sample}_virus_summary.tsv",
            sample=Samples.keys(),
        ),
        expand(
            VIRUS_FP / "cenote_taker" / "{sample}" / "{sample}_virus_sequences.fna",
            sample=Samples.keys(),
        ),
        VIRUS_FP / "cenote_taker" / "cenote_virus_summary.tsv",
        VIRUS_FP / "cenote_taker" / "cenote_virus_genomes.fna",
        VIRUS_FP / "cenote_taker" / "cenote_vOTU_repr_contig_map.tsv",
        VIRUS_FP / "cenote_taker" / "cenote_vOTU_repr_contigs.fna",
        VIRUS_FP / "cenote_taker" / "cenote_vOTU_mapped_read_counts.tsv",
        VIRUS_FP / "cenote_taker" / "cenote_vOTU_taxonomy.tsv",
        expand(
            VIRUS_FP / "cenote_taker" / "{sample}_reads_to_vOTUs.idxstats.tsv",
            sample=Samples.keys(),
        ),
        expand(
            VIRUS_FP / "cenote_taker" / "{sample}_reads_to_vOTUs.genomecov.tsv",
            sample=Samples.keys(),
        ),
#        expand(
#            VIRUS_FP / "cenote_taker" / "untranslated_blastx" / "{sample}.btf",
#            sample=Samples.keys(),
#        ),
#        expand(
#            VIRUS_FP / "cenote_taker" / "{sample}" / "ct_processing" / "final_taxonomy" / "virus_taxonomy_summary.tsv",
#            sample=Samples.keys(),
#        ),
#        expand(
#            VIRUS_FP / "cenote_taker" / "{sample}" / "ct_processing" / "mapping_reads" / "oriented_hallmark_contigs.pruned.coverage.tsv",
#            sample=Samples.keys(),
#        ),

rule run_cenote_taker:
    input:
        contigs=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        contigs=VIRUS_FP / "cenote_taker" / "{sample}" / "{sample}_virus_sequences.fna",
        ct3log=VIRUS_FP / "cenote_taker" / "{sample}" / "{sample}_cenotetaker.log",
        summary=VIRUS_FP / "cenote_taker" / "{sample}" / "{sample}_virus_summary.tsv",
    benchmark:
        BENCHMARK_FP / "cenote_taker_{sample}.tsv"
    log:
        LOG_FP / "cenote_taker_{sample}.log",
    params:
        out_dir=str(VIRUS_FP / "cenote_taker"),
        sample="{sample}",
        db_fp=Cfg["sbx_cenote_taker"]["cenote_taker_db"],
        molecule_type=Cfg["sbx_cenote_taker"]["molecule_type"],
        lin_minimum_hallmark_genes=Cfg["sbx_cenote_taker"]["lin_minimum_hallmark_genes"],
        seqtech=Cfg["sbx_cenote_taker"]["seqtech"],
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
        cenotetaker3 --contigs {input.contigs} --reads {input.r1} {input.r2} --seqtech {params.seqtech} --molecule_type {params.molecule_type} -r {params.sample} -p T --lin_minimum_hallmark_genes {params.lin_minimum_hallmark_genes} >> {log} 2>&1 

        if [ -s {output.ct3log} ] && [ ! -s {output.summary} ] && [ ! -s {output.contigs} ]; then
            if grep -q "no contigs with at least [0-9]\\+ hallmark genes of type(s) (virion rdrp) found" {output.ct3log}; then
                echo "Making dummy results for this sample because cenote-taker3 did not find enough hallmark genes in any contig:" >> {log}
                touch {output.contigs} {output.summary}
                echo "Exiting with exit code 0." >> {log}
                exit 0
            elif grep -q "couldn't find .*contigs_over_1000nt.fasta" {output.ct3log}; then
                echo "Making dummy results for this sample because cenote-taker3 did not find any viral contigs:" >> {log}
                touch {output.contigs} {output.summary}
                echo "Exiting with exit code 0." >> {log}
                exit 0
            fi
        elif [ -s {output.ct3log} ] && [ -s {output.summary} ] && [ -s {output.contigs} ]; then
            echo "Cenote-Taker3 log, summary, and contigs were generated." >> {log}
        else
            echo "Cenote-Taker3 log file is empty or does not exist" >> {log}
        fi
        """

rule ct3_output_merge:
    input:
        tables=expand(
            VIRUS_FP / "cenote_taker" / "{sample}" / "{sample}_virus_summary.tsv",
            sample=Samples.keys(),
        ),
        seqs=expand(
            VIRUS_FP / "cenote_taker" / "{sample}" / "{sample}_virus_sequences.fna",
            sample=Samples.keys(),
        ),
    output:
        tables=VIRUS_FP / "cenote_taker" / "cenote_virus_summary.tsv",
        seqs=VIRUS_FP / "cenote_taker" / "cenote_virus_genomes.fna",
    threads: 1
    log:
        LOG_FP / "ct3_output_merge.log",
    conda:
        "envs/cenote_taker_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-cenote-taker"
    shell:
        """ 
        awk 'FNR==1 {{ if (NR==1) print; next }} 1' {input.tables} > {output.tables}
        cat {input.seqs} > {output.seqs}
        """

rule ct3_cluster:
    input:
        VIRUS_FP / "cenote_taker" / "cenote_virus_genomes.fna"
    output:
        fltr=temp(VIRUS_FP / "cenote_taker" / "cenote_vOTU_fltr.txt"),
        ani=VIRUS_FP / "cenote_taker" / "cenote_virus_genomes_ani.tsv",
        aniids=VIRUS_FP / "cenote_taker" / "cenote_virus_genomes_ani.ids.tsv",
        map=VIRUS_FP / "cenote_taker" / "cenote_vOTU_repr_contig_map.tsv",
        size=VIRUS_FP / "cenote_taker" / "cenote_vOTU_num_contigs.tsv",
        seqs=VIRUS_FP / "cenote_taker" / "cenote_vOTU_repr_contigs.fna",
        bwaidx=VIRUS_FP / "cenote_taker" / "cenote_vOTU_repr_contigs.fna.amb",
    threads: 16
    resources:
        mem_mb=96000,
        runtime=1440,
        cpus_per_task=16,
    log:
        LOG_FP / "ct3_cluster_vOTUs.log",
    conda:
        "envs/vclust_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-cenote-taker"
    shell:
        """
        # Create a pre-alignment filter.
        vclust prefilter -i {input} -o {output.fltr} --min-ident 0.95 -t {threads} >> {log} 2>&1 
        # Calculate ANI measures for genome pairs specified in the filter.
        vclust align -i {input} -o {output.ani} --filter {output.fltr} -t {threads} >> {log} 2>&1 
        # Cluster contigs into vOTUs using the MIUVIG thresholds and the Leiden algorithm.
        # vclust cluster -i {output.ani} -o {output.map} --ids {output.aniids} --algorithm leiden --metric ani --ani 0.95 --qcov 0.85 >> {log} 2>&1 
        # Cluster for representative contigs
        vclust cluster -i {output.ani} -o {output.map} --out-repr --ids {output.aniids} --algorithm leiden --metric ani --ani 0.95 --qcov 0.85 >> {log} 2>&1 
        # Get only the representative clusters and count the number of contigs in the cluster
        echo -e "cluster\tcontig_count"; awk -F'\t' 'NR>1 {{count[$2]++}} END {{for (c in count) print c "\t" count[c]}}' {output.map} | sort -k2,2nr > {output.size}
        # Put all the sequences of representative clusters into one fasta file by extracting matching headers
        seqkit grep -f <(tail -n +2 {output.size} | cut -f1) {input} --threads {threads} > {output.seqs} 2>> {log}
        bwa index {output.seqs} >> {log} 2>&1 
        """

rule ct3_read_mapping:
    input:
        seqs=VIRUS_FP / "cenote_taker" / "cenote_vOTU_repr_contigs.fna",
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        bam=VIRUS_FP / "cenote_taker" / "{sample}_reads_to_vOTUs.bam",
        bai=VIRUS_FP / "cenote_taker" / "{sample}_reads_to_vOTUs.bam.bai",
        stats=VIRUS_FP / "cenote_taker" / "{sample}_reads_to_vOTUs.idxstats.tsv",
        cov=VIRUS_FP / "cenote_taker" / "{sample}_reads_to_vOTUs.genomecov.tsv",
    conda:
        "envs/vclust_env.yml"
    threads: 4
    resources:
        cpus_per_task=4,
    log:
        LOG_FP / "ct3_read_mapping_{sample}.log",
    shell:
        """
        bwa mem -t {threads} {input.seqs} {input.r1} {input.r2} | samtools view -u - | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam} {output.bai} >> {log} 2>&1 
        samtools idxstats {output.bam} | awk -v sampid={wildcards.sample} 'BEGIN{{print "ref\tseqlen\t" sampid "\tunmapped"}} {{print}}'> {output.stats} 2>> {log}
        bedtools genomecov -ibam {output.bam} -bga | awk 'BEGIN{{print "genome\tposition_start\tposition_end\tdepth"}} {{print}}'> {output.cov} 2>> {log}
        """

rule ct3_vOTU_table:
    input:
        stats=expand(
            VIRUS_FP / "cenote_taker" / "{sample}_reads_to_vOTUs.idxstats.tsv",
            sample=Samples.keys(),
        ),
        summary=VIRUS_FP / "cenote_taker" / "cenote_virus_summary.tsv",
    output:
        otutab=VIRUS_FP / "cenote_taker" / "cenote_vOTU_mapped_read_counts.tsv",
        taxtab=VIRUS_FP / "cenote_taker" / "cenote_vOTU_taxonomy.tsv",
    threads: 1
    log:
        LOG_FP / "ct3_vOTU_table.log",
    conda:
        "envs/cenote_taker_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-cenote-taker"
    shell:
        """ 
        # grab mapped reads and reference length for vOTU table
        paste {input.stats} | awk '{{printf "%s\t%s", $1, $2; for(i=3;i<=NF;i+=4) printf "\t%s", $i; print ""}}' > {output.otutab} 2>> {log}
        # grab the vOTU representative sequence taxonomic assignments by contig ID and put them in their own file
        #head -n1 {input.summary} > {output.taxtab}
        #grep -F "$(awk '{{print $1 "\\t"}}' {output.otutab})" {input.summary} 1>> {output.taxtab} 2>> {log}
        awk 'BEGIN {{FS=OFS="\\t"}} NR==FNR {{if (NR==1) {{h2=$0; next}}; right[$1]=$0; next}} FNR==1 {{print $0, h2; next}} {{if ($1 in right) {{print $0, right[$1]}} else {{print $0, ""}}}}' {input.summary} {output.otutab} > {output.taxtab} 2>> {log}
        """

# TODO: unclassified BLAST search against viral genome DB
#rule virus_blastx:
#    """Run blastx on untranslated genes against a target db and write to blast tabular format."""
#    input:
#        VIRUS_FP / "cenote_taker" / "{sample}" / "{sample}_virus_sequences.fna",
#    output:
#        VIRUS_FP / "cenote_taker" / "untranslated_blastx" / "{sample}.btf7.tsv",
#    benchmark:
#        BENCHMARK_FP / "run_virus_blastx_{sample}.tsv"
#    log:
#        LOG_FP / "run_virus_blastx_{sample}.log",
#    params:
#        blast_db=Cfg["sbx_cenote_taker"]["blast_db"],
#    threads: Cfg["sbx_cenote_taker"]["blastx_threads"]
#    resources:
#        mem_mb=24000,
#        runtime=720,
#    conda:
#        "envs/utils.yml"
#    container:
#        f"docker://sunbeamlabs/sbx_cenote_taker:{SBX_CENOTE_TAKER_VERSION}-utils"
#    shell:
#        """
#        if [ -s {input} ]; then
#            export BLASTDB=$(dirname {params.blast_db})
#            blastx \
#            -query {input} \
#            -db $(basename {params.blast_db}) \
#            -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
#            -num_threads {threads} \
#            -evalue 0.05 \
#            -max_target_seqs 100 \
#            -out {output} \
#            2>&1 | tee {log}
#        else
#            echo "Caught empty query" >> {log}
#            touch {output}
#        fi
#        """
