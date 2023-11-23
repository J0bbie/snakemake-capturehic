###############################################################################
# Rules relating to CHiCAGO analysis of Hi-C data.
###############################################################################

rule chicago_rmap:
    input:
        digest_file="genomes/{species}/digest_HindIII.txt",
    output:
        rmap="chicago/design/{species}/{species}.rmap",
    threads: 1
    resources:
        mem_mb=1024 * 2
    message: "Generating rmap for {wildcards.species}."
    shell:
        """
        tail -n +3 {input.digest_file} | awk 'BEGIN{{OFS=\"\t\";n=0}};{{n+=1; print $1, $2, $3, n }}' > {output.rmap}
        """

rule chicago_baitmap:
    input:
        rmap="chicago/design/{species}/{species}.rmap",
        probes=config["dir_supplementary"]+"probes/{species}_filt_finalTargets.bed"
    output:
        baitmap="chicago/design/{species}/{species}.baitmap",
        baitmap_4col="chicago/design/{species}/{species}.baitmap_4col.txt",
    conda:
        "envs/chicago.yaml",
    threads: 1
    resources:
        mem_mb=1024 * 2
    message: "Generating baitmap for {wildcards.species}."
    shell:
        """
        bedtools intersect -a {input.rmap}  -b {input.probes} -wa -wb | \
        bedtools groupby -g 1-4 -c 9 -o distinct > {output.baitmap}

        cat {output.baitmap} | cut -f1-4 > {output.baitmap_4col}
        """


rule chicago_designfiles:
    input:
        rmap="chicago/design/{species}/{species}.rmap",
        baitmap="chicago/design/{species}/{species}.baitmap",
        baitmap_4col="chicago/design/{species}/{species}.baitmap_4col.txt"
    output:
        f1="chicago/design/{species}/{species}.poe",
        f2="chicago/design/{species}/{species}.npb",
        f3="chicago/design/{species}/{species}.nbpb"
    conda:
        "envs/chicago.yaml",
    threads: 1
    resources:
        mem_mb=1024 * 5,
    params:
        designDir="chicago/design/{species}/"
    message: "Generating design files for {wildcards.species}."
    shell:
        """
        python2.7 $(which makeDesignFiles.py) \
                    --rmapfile={input.rmap}   \
                    --baitmapfile={input.baitmap} \
                    --outfilePrefix=chicago/design/{wildcards.species}/{wildcards.species} \
                    --minFragLen=150 \
                    --maxFragLen=40000 \
                    --designDir={params.designDir}
        """

rule chicago_output:
    input:
        bam="alignment/{species}/{sample}/{sample}_r_1_2.hicup.bam",
        rmap="chicago/design/{species}/{species}.rmap",
        baitmap="chicago/design/{species}/{species}.baitmap",
    output:
        chicago_input="chicago/input/{species}/{sample}.chinput",
        chicago_bedpe="chicago/input/{species}/{sample}_bait2bait.bedpe"
    params:
        out_name="chicago/input/{species}/{sample}"
    conda:
        "envs/chicago.yaml",
    message: "Generating CHiCAGO input files for ({wildcards.species}): {wildcards.sample}."
    shell:
        """

        mkdir -p chicago/input/{wildcards.species}/

        {workflow.basedir}/rules/misc/bam2chicago_v2.sh \
            --bamfile {input.bam} \
            --baitmap {input.baitmap} \
            --rmap {input.rmap} \
            --outname {params.out_name}
        
        # Move one folder back.
        mv {params.out_name}/* chicago/input/{wildcards.species}/
        rmdir {params.out_name}
        """

def get_input_tissue(wildcards):
    """"
    Get the replicates for a given tissue / species.
    """
    samples = samplesheet.query("species == @wildcards.species & tissue == @wildcards.tissue")["sample"].tolist()
    return samples
    
    
rule run_chicago:
    input:
        f1="chicago/design/{species}/{species}.poe",
        f2="chicago/design/{species}/{species}.npb",
        f3="chicago/design/{species}/{species}.nbpb",
        chicago_input=lambda w: expand("chicago/input/{species}/{sample}.chinput", species=w.species, sample=get_input_tissue(w)),
        chicago_bedpe=lambda w: expand("chicago/input/{species}/{sample}_bait2bait.bedpe", species=w.species, sample=get_input_tissue(w))
    output:
        chicago_ibed="chicago/output/{species}/{tissue}_chicago.ibed",
        chicago_seqmonk="chicago/output/{species}/{tissue}_chicago_seqmonk.txt",
        chicago_washU="chicago/output/{species}/{tissue}_chicago_washU_text.txt",
        chicago_rds="chicago/output/{species}/{tissue}_chicago.rds",
        chicago_full_table="chicago/output/{species}/{tissue}_chicago_full_table.txt",
    params:
        inputfolder="chicago/input/{species}/",
        designfolder="chicago/design/{species}/",
        outfolder="chicago/output/{species}/",
        samples=lambda w: get_input_tissue(w)
    conda:
        "envs/chicago.yaml",
    message: "Running CHiCAGO for {wildcards.species} x {wildcards.tissue}."
    shell:
        """
        {workflow.basedir}/rules/scripts/1.run_chicago.R \
        --inputfolder {params.inputfolder} \
        --designfolder {params.designfolder} \
        --samples {params.samples} \
        --species {wildcards.species} \
        --tissue {wildcards.tissue} \
        --out={params.outfolder}
        """
