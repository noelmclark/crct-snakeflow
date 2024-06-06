# these are rules that serve the purpose of running the
# workflow to generate particular files.  They are particularly
# here so that we can easily call for the generation of files
# that need to be made so that we can check things and set parameters
# for later in the workflow.  For example, we might want to know
# the distribution of QUAL and QD values so that we can set
# a cutoff for bootstrapping the BQSR.



# I need to fix these to be the actual paths for the given outputs  Will do that later...
#rule dest_qc:
#	input:
#		expand("results/qc/multiqc.html"),
##		expand("results/qc/bcftools_stats/all-pass-maf-{maf}.txt", maf=mafs),




# This turns chromosomes.tsv and scaffolds.tsv into scatter_intervals.tsv.
# Needs to have R installed and on the path.
rule dest_scatter_intervals:
    input:
        chroms=config["chromosomes"],
        scaffs=config["scaffold_groups"],
    params:
        binsize="{binsize}"
    envmodules: "R/4.2.2"
    output:
        tsv="results/scatter_config/scatters_{binsize}.tsv"
    log:
    	"results/logs/dest_scatter_intervals/log_{binsize}.txt"
    script:
    	"../scripts/sequence-scatter-bins.R"


# this is for downsampling the bams and nothing more.  If you want
# to downsample bams and then run through the entire gVCF workflow with those
# you should see the next rule...
#rule dest_downsample_bams_only:
#	input:
#		bam=expand(
#			"results/downsample-{cov}X/overlap_clipped/{sample}.bam",  
#			cov = config["downsample_bams"]["depths"],
#			sample = sample_list )



# this is just here to make it easy to do a run that just
# force-calls the sites
#rule force_call_sites:
#	input:
#		vcf="results/force-call/final.vcf.gz"
