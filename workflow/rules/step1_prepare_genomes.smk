# Rules relating to downloading the reference genomes and preparing them for subsequent analysis.
rule download_species:
    output:
        genome="genomes/{species}/genome.fa.gz",
    threads: 1
    resources:
        mem_mb=1024*5
    params:
        species=lambda wildcards: wildcards.species,
    conda:
        "envs/alignment.yaml"
    message: "Downloading reference genome of {wildcards.species}"
    shell:
        """
        if [ {params.species} == "canis_familiaris" ]; then
            for chr in {{1..38}} X MT; do
                wget -P genomes/{params.species}/ https://ftp.ensembl.org/pub/release-98/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna.chromosome.$chr.fa.gz;
            done
        fi

        if [ {params.species} == "macaca_mulatta" ]; then
            for chr in {{1..20}} X Y; do
                wget -P genomes/{params.species}/ https://ftp.ensembl.org/pub/release-98/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.$chr.fa.gz;
            done
        fi

        if [ {params.species} == "mus_musculus" ]; then
            for chr in {{1..19}} X Y MT; do
                wget -P genomes/{params.species}/ https://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.$chr.fa.gz;
            done
        fi

        if [ {params.species} == "rattus_norvegicus" ]; then
            for chr in {{1..20}} X Y MT; do
                wget -P genomes/{params.species}/ https://ftp.ensembl.org/pub/release-98/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.chromosome.$chr.fa.gz;
            done
        fi

        if [ {params.species} == "callithrix_jacchus" ]; then
            for chr in {{1..20}} X Y MT; do
                wget -P genomes/{params.species} https://ftp.ensembl.org/pub/release-98/fasta/callithrix_jacchus/dna/Callithrix_jacchus.ASM275486v1.dna.toplevel.fa.gz;
            done
        fi

        # Combine separate chromosomes into a single genome file.
        zcat genomes/{params.species}/*.fa.gz | seqkit sort --natural-order | gzip -c > genomes/{params.species}/genome.fa.gz

        # Remove separate chromosomes.
        find genomes/{params.species}/ -type f ! -name 'genome.fa.gz' -delete
        """

rule generate_digest:
    input:
        genome="genomes/{species}/genome.fa.gz",
    output:
        digest=directory("genomes/{species}/digest/"),
    threads: 20
    resources:
        mem_mb=1024*30
    conda:
        "envs/alignment.yaml"
    message: "Generating digest for {wildcards.species}"
    shell:
        """
        hicup_digester \
        --genome {wildcards.species} \
        -re1 A^AGCTT,HindIII \
        --outdir {output.digest} \
        {input.genome} 
        """


rule generate_bowtie2_index:
    input:
        genome="genomes/{species}/genome.fa.gz",
    output:
        idx=directory("genomes/{species}/bowtie2_index/"),
    threads: 20
    resources:
        mem_mb=1024*30
    conda:
        "envs/alignment.yaml"
    message: "Generating bowtie2 index for {wildcards.species}"
    shell:
        """
        mkdir -p genomes/{wildcards.species}/bowtie2_index
        bowtie2-build --threads {threads} {input.genome} {output.idx}/{wildcards.species}
        """
