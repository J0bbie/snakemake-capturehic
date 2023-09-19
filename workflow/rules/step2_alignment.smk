rule alignment_hicup:
    input:
        fastq = lambda w:
            [
                "{fastq}/{sample}_r_1.fq.gz".format(fastq=config['dir_input'], sample = w.sample),
                "{fastq}/{sample}_r_2.fq.gz".format(fastq=config['dir_input'], sample = w.sample),
            ],
        bowtie2_index = directory("genomes/{species}/bowtie2_index/"),
    output:
        out_dir = directory("alignment/{species}/{sample}/"),
        bam = "alignment/{species}/{sample}.bam",
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
        hicup  \
        --digest {params.digest} \
        --format Sanger \
        --index {input.bowtie2_index} \
        --longest 700 \ 
        --zip 1 \
        --keep 1 \
        --outdir  {output.out_dir} \
        --shortest 50 \
        --threads {threads} \
        {input.fastq[0]} {input.fastq[1]}
        """