### following rules are to run PSMC using cleaned BAM files
## based on lh3 documentation at: https://github.com/lh3/psmc
## I was curious if there would be a difference in estimates using the bcftools call (from-bam) vs gatk variant pipelines (from-bcf)

##rule to remove the y chrom (NC_048593.1) and mitogenome (NC_001717.1) and all scaffold regions (NW.*)
rule get_aut_regions:
    input:
        scat_path="results/scatter_config/scatters_1200000.tsv"
    output:
        regions="results/psmc/get-aut-regions/autosomal_regions.bed"
    log:
        "results/logs/psmc/get-aut-regions.log"
    benchmark:
        "results/benchmarks/psmc/get-aut-regions.bmk"
    shell:
        " (awk -v chrom='NC_048593.1' -f workflow/scripts/PSMC/remove_chrom_regions.awk {input.scat_path} | "
        " awk -v chrom='NC_001717.1' -f workflow/scripts/PSMC/remove_more_regions.awk | "
        " awk -v chrom='NW_.*' -f workflow/scripts/PSMC/remove_more_regions.awk "
        " > {output.regions} ) 2> {log} "

### run PSMC from BAMs
## rule to extract only regions that map to autosomes (not Y chrom or scaffolds)
# sex chroms are shown to have impacts on PSMC curves
# the y-chrom for the o. mykiss reference is NC_048593.1
rule get_aut_bams:
    input:
        regions="results/psmc/get-aut-regions/autosomal_regions.bed",
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",
        bai="results/angsd_bams/overlap_clipped/{sample}.bai"
    output:
        aut_bam="results/psmc/aut-bams/aut_{sample}.bam",
        aut_bai="results/psmc/aut-bams/aut_{sample}.bai"
    conda:
        "../envs/sambcftools.yaml"
    log:
        rem="results/logs/psmc/aut-bams/aut-{sample}.log",
        index="results/logs/psmc/aut-bams/index-{sample}.log",
    benchmark:
        "results/benchmarks/psmc/aut-bams/{sample}.bmk"
    shell:
        " samtools view -b -h -L {input.regions} -o {output.aut_bam} {input.bam} 2> {log.rem} && "
        " samtools index {output.aut_bam} -o {output.aut_bai} 2> {log.index} "


## rule to get a consensus fastq sequence file for PSMC
# uses the bam files with only reads mapping to autosomal chromosomes
# option -C 50 downgrades mapping quality (by coeff given) for reads containing excessive mismatches
# option -d sets and minimum read depth and -D sets the maximum 
# the flycatcher paper recommends setting -d 10 for depth coverage >10x
rule psmc_consensus_sequence:
    input:
        bam="results/psmc/aut-bams/aut_{sample}.bam",  
        ref="resources/genome/OmykA.fasta",
    output:
        "results/psmc/psmc-consensus-sequence/{sample}.fq.gz",
    conda:
        "../envs/sambcftools.yaml"
    resources:
        time="23:59:59",
        mem_mb=9400,
        cpus=2,
    threads: 2
    log:
        "results/logs/psmc/psmc-consensus-sequence/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/psmc-consensus-sequence/{sample}.bmk"
    shell:
        "bcftools mpileup --full-BAQ -C50 -Ou -f {input.ref} {input.bam} | bcftools call -c - | " 
        "vcfutils.pl vcf2fq -d 10 -D 36 | gzip > {output} 2> {log}"

# rule to create psmcfa file per sample
rule psmcfa:
    input:
        "results/psmc/psmc-consensus-sequence/{sample}.fq.gz"
    output:
        "results/psmc/psmcfa/{sample}.psmcfa"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/psmcfa/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/psmcfa/{sample}.bmk"
    shell:
        "fq2psmcfa -q20 {input} > {output} 2> {log}"

## rule to run psmc with approved options
rule run_psmc:
    input:
        "results/psmc/psmcfa/{sample}.psmcfa"
    output:
        "results/psmc/run-psmc/{sample}.psmc"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/run-psmc/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/run-psmc/{sample}.bmk"
    shell:
        "psmc -N25 -t10 -r5 -p '10+6*2+18*1+8*2+8*1' -o {output} {input} 2> {log}"





### PSMC bootstrapping ###
## these rules deal with running 100 bootstrap replicates of PSMC for each sample 

# this rule splits long chrom sequences to shorter segments for each sample
rule split_psmcfa:
    input:
        "results/psmc/psmcfa/{sample}.psmcfa"
    output:
        "results/psmc/bootstrap/split-psmcfa/{sample}-split.psmcfa"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/bootstrap/split-psmcfa/{sample}-split.psmcfa.log"
    benchmark:
        "results/benchmarks/psmc/bootstrap/split-psmcfa/{sample}-split.psmcfa.bmk"
    shell:
        " splitfa {input} > {output} 2> {log} "

# this rule runs PSMC 100 times with a random sample (-b) of short chroms from the split file each of 100 times (wildcard.psmcbootround=1:100)
rule bootstrap_psmc:
    input:
        split="results/psmc/bootstrap/split-psmcfa/{sample}-split.psmcfa"
    output:
        "results/psmc/bootstrap/run-psmc/{sample}/{sample}-bootround-{psmcbootround}.psmc"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/bootstrap/run-psmc/{sample}/{sample}-bootround-{psmcbootround}.log"
    benchmark:
        "results/benchmarks/psmc/bootstrap/run-psmc/{sample}/{sample}-bootround-{psmcbootround}.bmk"
    shell:
        "psmc -N25 -t10 -r5 -b -p '10+6*2+18*1+8*2+8*1' -o {output} {input.split} 2> {log} "

# this rule merges the full PSMC and all PSMC results of 100 bootstraps for each sample
rule merge_psmc_boostraps:
    input:
        whole="results/psmc/run-psmc/{sample}.psmc",
        boots=expand("results/psmc/bootstrap/run-psmc/{{sample}}/{{sample}}-bootround-{pbr}.psmc", pbr=psmcbootround)
    output:
        "results/psmc/bootstrap/run-psmc/{sample}/{sample}-100bootstrap.psmc"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/bootstrap/run-psmc/{sample}/{sample}-100bootstrap.log"
    benchmark:
        "results/benchmarks/psmc/bootstrap/run-psmc/{sample}/{sample}-100bootstrap.bmk"
    shell:
        " cat {input.whole} > {output} && " #put whole psmc file into output file first
        " for i in {input.boots}; do "
        " (cat $i); " # then add each of the psmc bootstrap outputs
        " done >> {output} 2> {log} "





### PSMC plotting ###

## rule to plot psmc for each individual
# apparently, you can set it to multiline mode using -M and supplying all psmc files you want plotted together
# can also set min (-x) and max (-X) generations for mapping 
# -u [per-generation mutation rate] from https://doi.org/10.1371/journal.pgen.1010918
# -g [generation time in years] 
rule psmc_plot:
    input:
        "results/psmc/run-psmc/{sample}.psmc"
    output:
        eps="results/psmc/psmc-plot/by-sample/{sample}.eps",
        par="results/psmc/psmc-plot/by-sample/{sample}.par"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/psmc-plot/by-sample/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/psmc-plot/by-sample/{sample}.bmk"
    shell:
        "psmc_plot.pl -u 8.0e-09 -g 3 {output.eps} {input} 2> {log}"


## rule to plot all from-bam PSMC by subsamp code - use this one the most
# I wrote this up to easily get overlain plots with different combinations of individuals 
rule psmc_plot_by_subsamp:
    input:
        psmc=get_psmc_subsamps,
    params:
        samps=get_comma_sep_subsamp_names,
    output:
        "results/psmc/psmc-plot/{subsamp}/{subsamp}",
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/psmc-plot/{subsamp}.log"
    benchmark:
        "results/benchmarks/psmc/psmc-plot/{subsamp}.bmk"
    shell:
        " psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" -M {params.samps} {output} {input.psmc} 2> {log} "


rule psmc_plot_by_subsamp_upper_bound:
    input:
        psmc=get_psmc_subsamps,
    params:
        samps=get_comma_sep_subsamp_names,
    output:
        "results/psmc/psmc-plot/{subsamp}-upper-bound/{subsamp}-upper-bound",
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/psmc-plot/{subsamp}-upper-bound.log"
    benchmark:
        "results/benchmarks/psmc/psmc-plot/{subsamp}-upper-bound.bmk"
    shell:
        " psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" -Y 10 -M {params.samps} {output} {input.psmc} 2> {log} "



## rule to plot all PSMC outputs together! super messy so don't normally use
# can also set min (-x) and max (-X) generations for mapping
# -P sets the legend position based on gnuplot syntax: https://gnuplot.sourceforge.net/docs_4.2/node192.html 
# explanation of psmc_plot.pl options https://github.com/lh3/psmc/blob/master/utils/psmc_plot.pl
rule psmc_plot_all:
    input:
        psmc=expand("results/psmc/run-psmc/{s}.psmc", s=sample_list),
    params:
        samps=get_comma_sep_sample_names,
    output:
        "results/psmc/psmc-plot/all/all-together",
        #par="results/psmc/psmc-plot/all/all-together.par"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/psmc-plot/all/all-together.log"
    benchmark:
        "results/benchmarks/psmc/psmc-plot/all/all-together.bmk"
    shell:
        " psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" -M {params.samps} {output} {input.psmc} 2> {log}"


## rule to plot whole PSMC and PSMC 100 bootstraps per sample, yay
rule psmc_plot_bootstrap:
    input:
        "results/psmc/bootstrap/run-psmc/{sample}/{sample}-100bootstrap.psmc",
    output:
        "results/psmc/bootstrap/psmc-plot/{sample}-100bootstrap",
        #par="results/psmc/bootstrap/psmc-plot/{sample}-100bootstrap.par"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/bootstrap/psmc-plot/{sample}-100bootstrap.log"
    benchmark:
        "results/benchmarks/psmc/bootstrap/psmc-plot/{sample}-100bootstrap.bmk"
    shell:
        " psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" {output} {input} 2> {log} "

## rules to plot PSMC with 100 bootstraps per sample by subsamp -- not sure how it will look or if it will be useful
rule psmc_plot_bootstrap_by_subsamp:
    input:
        psmc=get_psmc_boot_subsamps,
    params:
        samps=get_comma_sep_subsamp_names,
    output:
        "results/psmc/bootstrap/psmc-plot/{subsamp}/{subsamp}",
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/bootstrap/psmc-plot/{subsamp}.log"
    benchmark:
        "results/benchmarks/psmc/bootstrap/psmc-plot/{subsamp}.bmk"
    shell:
        " psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" -M {params.samps} {output} {input.psmc} 2> {log} "





### PSMC Parameter Testing ###
## Here is the shell I used to test 70 different PSMC paramter patterns
## I varied -t and -p until I decided on a -t10 and -p '10+6*2+18*1+8*2+8*1'
# That -p pattern gave me at least 10 recombinations in the 20th round of psmc 
# in each atomic time interval (5th column of .psmc output file)
# Commented out so I don't have to save the iteration files when I'm done
# I kept track of patterns tested and outcomes in an external spreadsheet

## rule to test psmc parameters many options 
#rule run_psmc_param_test1:
#    input:
#        "results/psmc/psmcfa/C106394.psmcfa"
#    output:
#        "results/psmc/run-psmc/param-test/i68_C106394.psmc"
#    conda:
#        "../envs/psmc.yaml"
#    log:
#        "results/logs/psmc/run-psmc/param-test/i68_C106394.log"
#    benchmark:
#        "results/benchmarks/psmc/run-psmc/param-test/i68_C106394.bmk"
#    shell:
#        "psmc -N25 -t10 -r5 -p '10+6*3+30*1+8*2' -o {output} {input} 2> {log}"

#rule run_psmc_param_test2:
#    input:
#        "results/psmc/psmcfa/C106394.psmcfa"
#    output:
#        "results/psmc/run-psmc/param-test/i69_C106394.psmc"
#    conda:
#        "../envs/psmc.yaml"
#    log:
#        "results/logs/psmc/run-psmc/param-test/i69_C106394.log"
#    benchmark:
#        "results/benchmarks/psmc/run-psmc/param-test/i69_C106394.bmk"
#    shell:
#        "psmc -N25 -t10 -r5 -p '12+6*4+20*1+8*2' -o {output} {input} 2> {log}"

#rule run_psmc_param_test3:
#    input:
#        "results/psmc/psmcfa/C106394.psmcfa"
#    output:
#        "results/psmc/run-psmc/param-test/i70_C106394.psmc"
#    conda:
#        "../envs/psmc.yaml"
#    log:
#        "results/logs/psmc/run-psmc/param-test/i70_C106394.log"
#    benchmark:
#        "results/benchmarks/psmc/run-psmc/param-test/i70_C106394.bmk"
#    shell:
#        "psmc -N25 -t10 -r5 -p '10+6*2+18*1+8*2+8*1' -o {output} {input} 2> {log}"

## This rule plots the different patterns for one sample together to compare 
#rule psmc_plot_params:
#    input:
#        psmc65="results/psmc/run-psmc/param-test/i68_C106394.psmc",
#        psmc69="results/psmc/run-psmc/param-test/i69_C106394.psmc",
#        psmc70="results/psmc/run-psmc/param-test/i70_C106394.psmc",
#    output:
#        "results/psmc/psmc-plot/param-test/param-test-2",
#        #par="results/psmc/psmc-plot/param-test/param-test.par"
#    conda:
#        "../envs/psmc.yaml"
#    log:
#        "results/logs/psmc/psmc-plot/param-test/param-test-2.log"
#    benchmark:
#        "results/benchmarks/psmc/psmc-plot/param-test/param-test-2.bmk"
#    shell:
#        " psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" -M i68,i69,i70 {output} {input.psmc68} {input.psmc69} {input.psmc70} 2> {log}"