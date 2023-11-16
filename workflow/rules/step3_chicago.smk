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
        baitmap_4col="chicago/design/{species}.baitmap_4col.txt",
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
        rmap="chicago/design/{species}.rmap",
        baitmap="chicago/design/{species}.baitmap",
        baitmap_4col="chicago/design/{species}.baitmap_4col.txt"
    output:
        f1="chicago/design/{species}.poe",
        f2="chicago/design/{species}.npb",
        f3="chicago/design/{species}.nbpb"
    conda:
        "envs/chicago.yaml",
    message: "Generating design files for {wildcards.species}."
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

# # Combine sample replicates into a single file.
# def make_table(sample_list, otp_Dir, species):
#     samples_of_species = list(filter(lambda x: species in x, sample_list))
#     flname = otp_Dir + "peakMatrix_f/" + species + "_toPeakMatrix.tab"

#     for samp in samples_of_species:
#         sp = samp.split("_")[2]
#         tiss = samp.split("_")[3]
#         combo = sp + "_" + tiss
#         row = samp + "\t" + otp_Dir + combo + "/" + samp + ".rds"
       
#         with open(flname, "a") as otp:
#             otp.write(row + "\n")
#     return

# rule create_matrixInput:
# 	output:
# 		otp_fl = "chicago/input/{species}_toPeakMatrix.tab"
# 	params:
# 		sample_list = samplesheet.query("species == '{species}'")["sample"].tolist(),
# 	run:
#         make_table(sample_list = params.sample_list, otp_Dir = "chicago/design/", species = wildcards.species)
