###############################################################################
# Rules relating to CHiCAGO analysis of Hi-C data.
###############################################################################

rule chicago_rmap:
    input:
        digest_file="genomes/{species}/digest_HindIII.txt",
    output:
        rmap="chicago/design/{species}.rmap",
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
        rmap="chicago/design/{species}.rmap",
        probes=config["dir_supplementary"]+"probes/{species}_filt_finalTargets.bed"
    output:
        baitmap="chicago/design/{species}.baitmap",
    conda:
        "envs/chicago.yaml",
    threads: 1
    resources:
        mem_mb=1024 * 2
    shell:
        """
        bedtools intersect -a {input.rmap}  -b {input.probes} | \
        bedtools groupby -g 1-4 -c 4 -o distinct > {output.baitmap}
        """

rule chicago_designfiles:
    input:
        rmap="chicago/design/{species}.rmap",
        baitmap="chicago/design/{species}.baitmap",
    output:
        f1="chicago/design/{species}.poe",
        f2="chicago/design/{species}.npb",
        f3="chicago/design/{species}.nbpb"
    conda:
        "envs/chicago.yaml",
    shell:
        """
        python2.7 $(which makeDesignFiles.py) \
                    --rmapfile={input.rmap}   \
                    --baitmapfile={input.baitmap} \
                    --outfilePrefix={wildcards.species} \
                    --minFragLen=150 \
                    --maxFragLen=40000 \
                    --designDir=chicago/design/
        """

rule chicago_output:
    input:
        bam="alignment/{species}/{sample}/{sample}_r_1_2.hicup.bam",
        rmap="chicago/design/{species}.rmap",
        baitmap="chicago/design/{species}.baitmap",
    output:
        chicago_input="chicago/design/{species}/{sample}.chicinput",
        chicago_bedpe="chicago/design/{species}/{sample}.bait2bait.bedpe"
    params:
        out_name="chicago/design/{species}/{sample}"
    conda:
        "envs/chicago.yaml",
    shell:
        """
        bam2chicago.sh \
            --bamfile {input.bam} \
            --baitmap {input.baitmap} \
            --rmap {input.rmap} \
            --outname {params.out_name}
        """
