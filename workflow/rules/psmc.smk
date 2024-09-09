### following rules are to run PSMC
## based on lh3 documentation at: https://github.com/lh3/psmc

##rule to test the remove chrom regions awk script
rule remove_y_regions:
    input:
        scat_path="results/scatter_config/scatters_1200000.tsv"
    output:
        regions="results/psmc/remove-y-regions/autosomal_regions.bed"
    log:
        "results/logs/psmc/remove-y-regions.log"
    benchmark:
        "results/benchmarks/psmc/remove-y-regions.bmk"
    shell:
        " awk -v chrom='NC_048593.1' -f workflow/scripts/PSMC/remove_chrom_regions.awk "
        " {input.scat_path} > {output.regions} 2> {log} "

## rule to split bam files into aut and y-chrom bam files
# sex chroms are shown to have impacts on PSMC curves
# the y-chrom for the o. mykiss reference is NC_048593.1
rule remove_sex_bams:
    input:
        regions="results/psmc/remove-y-regions/autosomal_regions.bed",
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",
        bai="results/angsd_bams/overlap_clipped/{sample}.bai"
    output:
        aut_bam="results/psmc/remove-sex-bams/aut_{sample}.bam",
        aut_bai="results/psmc/remove-sex-bams/aut_{sample}.bai"
    conda:
        "../envs/sambcftools.yaml"
    log:
        rem="results/logs/psmc/remove-sex-bams/remove-{sample}.log",
        index="results/logs/psmc/remove-sex-bams/index-{sample}.log",
    benchmark:
        "results/benchmarks/psmc/remove-sex-bams/{sample}.bmk"
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
        bam="results/psmc/remove-sex-bams/aut_{sample}.bam",  
        ref="resources/genome/OmykA.fasta",
    output:
        temp("results/psmc/psmc-consensus-sequence/{sample}.fq.gz") #temp removes files once they're no longer needed downstream
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
#rule psmcfa:
#    input:
#        "results/psmc/psmc-consensus-sequence/{sample}.fq.gz"
#    output:
#        "results/psmc/psmcfa/{sample}.psmcfa"
#    conda:
#        "../envs/psmc.yaml"
#    log:
#        "results/logs/psmc/psmcfa/{sample}.log"
#    benchmark:
#        "results/benchmarks/psmc/psmcfa/{sample}.bmk"
#    shell:
#        "fq2psmcfa -q20 {input} > {output} 2> {log}"






### PSMC plotting ###
## following rules are to test the effect of different parameter values on running PSMC
## namely the -t -p and -r in run_psmc

## rule to run psmc with whitefish options
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
        "psmc -N25 -t5 -r5 -p '4+20*2+6*4+4' -o {output} {input} 2> {log}"
        # could try splitting up some of the time intervals


## rule to plot psmc to visualize result
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


## rule to plot all PSMC outputs together!
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


## rule to plot all PSMC by subsamp code
# I wrote this up to easily get overlain plots with different combinations of individuals 
rule psmc_plot_by_subsamp:
    input:
        psmc=get_psmc_subsamps,
    params:
        samps=get_comma_sep_subsamp_names,
    output:
        "results/psmc/psmc-plot/by-{psmc_id}/{subsamp}/{subsamp}",
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/psmc-plot/by-{psmc_id}/{subsamp}.log"
    benchmark:
        "results/benchmarks/psmc/psmc-plot/by-{psmc_id}/{subsamp}.bmk"
    shell:
        " psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" -M {params.samps} {output} {input.psmc} 2> {log} "



### I was curious if there would be a difference in estimates using the bcftools call vs gatk variant pipelines
## so this chunk of rules runs PSMC from the BCF files

# first we need to break the BCF into samples and remove the variants mapping to the y chrom
rule psmc_consensus_seq_from_vcf:
    input:
        bcf="results/bcf/all.bcf",
        tbi="results/bcf/all.bcf.csi",
        regions="results/psmc/remove-y-regions/autosomal_regions.bed",
    output:
        "results/psmc/from-vcf/psmc-consensus-sequence/{sample}.fq.gz"
    conda:
        "../envs/sambcftools.yaml"
    log:
        "results/logs/psmc/from-vcf/psmc-consensus-sequence/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/from-vcf/psmc-consensus-sequence/{sample}.bmk"
    shell:
        " bcftools query -s {wildcards.sample} -R {input.regions} {input.bcf} | "
        " vcfutils.pl vcf2fq -d 10 -D 36 | gzip > {output} 2> {log} "






## rule to test psmc parameters many options 
rule run_psmc_param_test1:
    input:
        "results/psmc/psmcfa/C106394.psmcfa"
    output:
        "results/psmc/run-psmc/param-test/i21_C106394.psmc"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/run-psmc/param-test/i21_C106394.log"
    benchmark:
        "results/benchmarks/psmc/run-psmc/param-test/i21_C106394.bmk"
    shell:
        "psmc -N25 -t10 -r5 -p '4*1+13*2+4+6*5' -o {output} {input} 2> {log}"

rule run_psmc_param_test2:
    input:
        "results/psmc/psmcfa/C106394.psmcfa"
    output:
        "results/psmc/run-psmc/param-test/i22_C106394.psmc"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/run-psmc/param-test/i22_C106394.log"
    benchmark:
        "results/benchmarks/psmc/run-psmc/param-test/i22_C106394.bmk"
    shell:
        "psmc -N25 -t10 -r5 -p '4*1+10*2+40' -o {output} {input} 2> {log}"

rule run_psmc_param_test3:
    input:
        "results/psmc/psmcfa/C106394.psmcfa"
    output:
        "results/psmc/run-psmc/param-test/i23_C106394.psmc"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/run-psmc/param-test/i23_C106394.log"
    benchmark:
        "results/benchmarks/psmc/run-psmc/param-test/i23_C106394.bmk"
    shell:
        "psmc -N25 -t10 -r5 -p '4*1+6*10' -o {output} {input} 2> {log}"
        # could try splitting up some of the time intervals

rule run_psmc_param_test4:
    input:
        "results/psmc/psmcfa/C106394.psmcfa"
    output:
        "results/psmc/run-psmc/param-test/i24_C106394.psmc"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/run-psmc/param-test/i24_C106394.log"
    benchmark:
        "results/benchmarks/psmc/run-psmc/param-test/i24_C106394.bmk"
    shell:
        "psmc -N25 -t10 -r5 -p '4*1+2*13+1*10+3*8' -o {output} {input} 2> {log}"
        # could try splitting up some of the time intervals

rule run_psmc_param_test5:
    input:
        "results/psmc/psmcfa/C106394.psmcfa"
    output:
        "results/psmc/run-psmc/param-test/i25_C106394.psmc"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/run-psmc/param-test/i25_C106394.log"
    benchmark:
        "results/benchmarks/psmc/run-psmc/param-test/i25_C106394.bmk"
    shell:
        "psmc -N25 -t10 -r5 -p '4*1+1*60' -o {output} {input} 2> {log}"
        # could try splitting up some of the time intervals

rule run_psmc_param_test6:
    input:
        "results/psmc/psmcfa/C106394.psmcfa"
    output:
        "results/psmc/run-psmc/param-test/i26_C106394.psmc"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/run-psmc/param-test/i26_C106394.log"
    benchmark:
        "results/benchmarks/psmc/run-psmc/param-test/i26_C106394.bmk"
    shell:
        "psmc -N25 -t10 -r5 -p '4*1+1*60' -o {output} {input} 2> {log}"
        # could try splitting up some of the time intervals