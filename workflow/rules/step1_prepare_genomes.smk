###############################################################################
# Rules relating downloading and processing the reference genomes.
###############################################################################

rule download_species:
    output:
        genome="genomes/{species}/genome.fa.gz",
    threads: 1
    resources:
        mem_mb=1024*20
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
            wget -P genomes/{params.species} https://ftp.ensembl.org/pub/release-98/fasta/callithrix_jacchus/dna/Callithrix_jacchus.ASM275486v1.dna.toplevel.fa.gz;

            # Rename and subset the chromosomes.
            cut -f1 {workflow.basedir}/rules/misc/callithrix_jacchus_scaffolds.txt > genomes/patterns.txt

            seqkit grep -n -f genomes/patterns.txt genomes/{params.species}/Callithrix_jacchus.ASM275486v1.dna.toplevel.fa.gz | gzip -c > genomes/{params.species}/genome_subsetted.fa.gz
            rm genomes/{params.species}/Callithrix_jacchus.ASM275486v1.dna.toplevel.fa.gz
            rm genomes/patterns.txt
        fi

        # Combine separate chromosomes into a single genome file.
        zcat genomes/{params.species}/*.fa.gz | seqkit sort --natural-order > genomes/{params.species}/genome.fa

        # Remove separate chromosomes.
        find genomes/{params.species}/ -type f ! -name 'genome.fa' -delete
        
        # Compress the genome file.
        pbgzip genomes/{params.species}/genome.fa
        """

rule download_features:
    output:
        gtf="genomes/{species}/{species}.gtf.gz",
    threads: 1
    params:
        species=lambda wildcards: wildcards.species,
    conda:
        "envs/alignment.yaml"
    message: "Downloading GTF of {wildcards.species}"
    shell:
        """
        if [ {params.species} == "canis_familiaris" ]; then
            wget -O {output} ftp://ftp.ensembl.org/pub/release-98/gtf/canis_familiaris/Canis_familiaris.CanFam3.1.98.gtf.gz;
        fi

        if [ {params.species} == "macaca_mulatta" ]; then
            wget -O {output} ftp://ftp.ensembl.org/pub/release-98/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.98.gtf.gz;
        fi

        if [ {params.species} == "mus_musculus" ]; then
            wget -O {output} ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz;
        fi

        if [ {params.species} == "rattus_norvegicus" ]; then
            wget -O {output} ftp://ftp.ensembl.org/pub/release-98/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.98.gtf.gz;
        fi
    
        if [ {params.species} == "callithrix_jacchus" ]; then
            wget -O {output} ftp://ftp.ensembl.org/pub/release-98/gtf/callithrix_jacchus/Callithrix_jacchus.ASM275486v1.98.gtf.gz;
        fi
        """

rule generate_chromsizes:
    input:
        genome="genomes/{species}/genome.fa.gz",
    output:
        fai="genomes/{species}/genome.fa.gz.fai",
        gzi="genomes/{species}/genome.fa.gz.gzi",
        chromsizes="genomes/{species}/genome.chrom.sizes",
    threads: 1
    resources:
        mem_mb=1024*10
    conda:
        "envs/alignment.yaml"
    message: "Calculating chromosome sizes for {wildcards.species}"
    shell:
        """
        samtools faidx {input.genome}
        cut -f1,2 {input.genome}.fai > {output.chromsizes}
        """


rule generate_annotations:
    input:
        gtf="genomes/{species}/{species}.gtf.gz",
        chromsizes="genomes/{species}/genome.chrom.sizes",
    output:
        bed_exons="genomes/{species}/{species}_exons.bed",
        bed_exons_merged="genomes/{species}/{species}_exons_merged.bed",
        bed_genes="genomes/{species}/{species}_genes.bed",
        bed_introns="genomes/{species}/{species}_introns.bed",
        bed_intergenic="genomes/{species}/{species}_intergenic.bed"
    threads: 1
    resources:
        mem_mb=1024*10
    conda:
        "envs/chicago.yaml"
    message: "Generating annotations for {wildcards.species}"
    shell:
        """
        chromosomes=$(cut -f1 {input.chromsizes})
        chromosomes2=\"($(echo ^$chromosomes | sed \'s/ /|^/g\'))\"

        # Extract exons.
        zcat {input.gtf} | awk 'BEGIN{{OFS="\\t"}};$3 == "exon" {{print $1,$4,$5}}' | grep -E $chromosomes2 | bedtools sort -g {input.chromsizes} -i - > {output.bed_exons}

        # Extract genes.
        zcat {input.gtf} | awk 'BEGIN{{OFS="\\t"}};$3 == "gene" {{print $1,$4,$5}}' | grep -E $chromosomes2 | bedtools sort -g {input.chromsizes} -i - > {output.bed_genes}

        # Extract introns.
        bedtools merge -i {output.bed_exons} > {output.bed_exons_merged}
        bedtools subtract -a {output.bed_genes} -b {output.bed_exons_merged} > {output.bed_introns}

        # Extract intergenic regions.
        bedtools complement -i {output.bed_genes} -g {input.chromsizes} > {output.bed_intergenic}
        """


rule generate_digest:
    input:
        genome="genomes/{species}/genome.fa.gz",
    output:
        digest_file="genomes/{species}/digest_HindIII.txt",
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
        --outdir genomes/{wildcards.species}/ \
        {input.genome} 

        # Rename the digest file.
        mv genomes/{wildcards.species}/Digest_*.txt {output.digest_file}
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
        mkdir -p {output.idx}
        bowtie2-build --threads {threads} {input.genome} {output.idx}/{wildcards.species}
        """
