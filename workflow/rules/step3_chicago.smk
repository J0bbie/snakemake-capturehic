###############################################################################
# Rules relating to CHiCAGO analysis of Hi-C data.
###############################################################################

rule chicago_rmap:
	input:
		"genomes/{species}/digest_HindIII.txt",
	output:
        "chicago/design/{species}.rmap",
    log:
        "logs/chicago/{species}_chicago_rmap.log",
    threads: 1
    resources:
        mem_mb=1024 * 2
	shell:
		"""
		mkdir -p chicago/design/{wildcards.species}
		tail -n +3 {input} | awk 'BEGIN{OFS=\"\t\";n=0};{n+=1; print \$1, \$2, \$3, n }' > {output}
		"""

rule chicago_baitmap:
	input:
		rmap = "chicago/design/{species}.rmap",
		probes = config["dir_supplementary"] + "{species}_filt_finalTargets.bed"
	output:
        "chicago/design/{species}.baitmap",
	conda:
		"envs/chicago.yaml",
	log:
		"logs/chicago/{species}_chicago_baitmap.log",
	threads: 1
	resources:
		mem_mb=1024 * 2
	shell:
		"""
		bedtools intersect -a {input.rmap}  -b {input.probes} | \
        bedtools groupby -i stdin -g 1-4 -c 9 -o distinct > {output}
		"""