# Rules relating to downloading the reference genomes and preparing them for subsequent analysis.
rule download_species:
    output:
        genome="genomes/{species}/genome.fa.gz",
    threads: 1
    resources:
        mem_mb=1024*5
    conda:
        "envs/alignment.yaml"
    message: "Downloading reference genome of {wildcards.species}"
    shell:
        """
        if {wildcards.species} == "canis_familiaris":
            for chr in {1..38} X MT; do
                wget -P genomes/{wildcards.species}/ https://ftp.ensembl.org/pub/release-98/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna.chromosome.$chr.fa.gz
            done

        if {wildcards.species} == "macaca_mulatta":
            for chr in {1..20} X Y; do
                wget -P genomes/{wildcards.species}/ https://ftp.ensembl.org/pub/release-98/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.$chr.fa.gz
            done
        if {wildcards.species} == "mus_musculus":
            for chr in {1..19} X Y MT; do
                wget -P genomes/{wildcards.species} https://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.$chr.fa.gz
            done

        if {wildcards.species} == "rattus_norvegicus":
            for chr in {1..20} X Y MT; do
                wget -P genomes/{wildcards.species} https://ftp.ensembl.org/pub/release-98/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.chromosome.$chr.fa.gz
            done

        if {wildcards.species} == "callithrix_jacchus":
            wget -P genomes/{wildcards.species} https://ftp.ensembl.org/pub/release-98/fasta/callithrix_jacchus/dna/Callithrix_jacchus.ASM275486v1.dna.toplevel.fa.gz

        # Combine all chromosomes
        zcat genomes/{wildcards.species}/*.fa.gz | seqkit sort --natural-order | gzip -c > genomes/{wildcards.species}/genome.fa.gz

        # Remove separate chromosomes.
        find genomes/{wildcards.species}/ -type f ! -name 'genome.fa.gz' -delete
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