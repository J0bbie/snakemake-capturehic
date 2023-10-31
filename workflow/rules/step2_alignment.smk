###############################################################################
# Rules relating to the alignment of Hi-C reads.
###############################################################################

rule alignment_hicup:
    input:
        fastq=lambda w:
            [
                "{fastq}{sample}_r_1.fq.gz".format(fastq=config['dir_input'], sample = w.sample),
                "{fastq}{sample}_r_2.fq.gz".format(fastq=config['dir_input'], sample = w.sample),
            ],
        genome_index="genomes/{species}/bowtie2_index/",
        genome_digest="genomes/{species}/digest_HindIII.txt",
    output:
        out_dir=directory("alignment/{species}/{sample}/"),
        bam="alignment/{species}/{sample}/{sample}_r_1_2.hicup.bam",
    log:
        "logs/alignment/{species}_{sample}_alignment_hicup.log",
    threads: 12
    resources:
        mem_mb=1024 * 25
    conda:
        "envs/alignment.yaml"
    message: "Running HiCUP for {wildcards.sample}."
    shell:
        """
        mkdir -p {output.out_dir}

        hicup  \
        --bowtie2 $(which bowtie2) \
        --digest {input.genome_digest} \
        --format Sanger \
        --index {input.genome_index}/{wildcards.species} \
        --longest 700 \
        --zip \
        --outdir {output.out_dir} \
        --shortest 50 \
        --threads {threads} \
        {input.fastq[0]} {input.fastq[1]} >& {log}
        """