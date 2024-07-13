## following rules are to test the effect of different parameter values on running PSMC
# namely the -d 10 in rule psmc_consensus_sequence 
# and the -t -p and -r in run_psmc
# based on lh3 documentation at: https://github.com/lh3/psmc

## rule to get a consensus fastq sequence file for PSMC
# option -C 50 downgrades mapping quality (by coeff given) for reads containing excessive mismatches
# option -d sets and minimum read depth and -D sets the maximum 
# the flycatcher paper recommends setting -d 10 for depth coverage >10x
rule psmc_consensus_sequence_test:
    input:
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",  
        ref="resources/genome/OmykA.fasta",
    output:
        "results/psmc-test/psmc-consensus-sequence-test/{sample}.fq.gz"
    conda:
        "../envs/sambcftools.yaml"
    resources:
        time="23:59:59",
        mem_mb=9400,
        cpus=2,
    threads: 2
    log:
        "results/logs/psmc-test/psmc-consensus-sequence-test/{sample}.log"
    benchmark:
        "results/benchmarks/psmc-test/psmc-consensus-sequence-test/{sample}.bmk"
    shell:
        "bcftools mpileup --full-BAQ -C50 -Ou -f {input.ref} {input.bam} | bcftools call -c - | " 
        "vcfutils.pl vcf2fq -d 10 -D 36 | gzip > {output} 2> {log}"


# rule to create psmcfa file per sample
rule psmcfa_test:
    input:
        "results/psmc-test/psmc-consensus-sequence-test/{sample}.fq.gz"
    output:
        "results/psmc-test/psmcfa-test/{sample}.psmcfa"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc-test/psmcfa-test/{sample}.log"
    benchmark:
        "results/benchmarks/psmc-test/psmcfa-test/{sample}.bmk"
    shell:
        "fq2psmcfa -q20 {input} > {output} 2> {log}"


# rule to run psmc
rule run_psmc_test:
    input:
        "results/psmc-test/psmcfa-test/{sample}.psmcfa"
    output:
        "results/psmc-test/run-psmc-test/{sample}.psmc"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc-test/run-psmc-test/{sample}.log"
    benchmark:
        "results/benchmarks/psmc-test/run-psmc-test/{sample}.bmk"
    shell:
        "psmc -N25 -t15 -r5 -p '4+25*2+4+6' -o {output} {input} 2> {log}"


## rule to plot psmc to visualize result
# apparently, you can set it to multiline mode using -M and supplying all psmc files you want plotted together
# can also set min (-x) and max (-X) generations for mapping 
# -u [per-generation mutation rate] from https://doi.org/10.1371/journal.pgen.1010918
# -g [generation time in years] 
rule psmc_plot_test:
    input:
        "results/psmc-test/run-psmc-test/{sample}.psmc"
    output:
        eps="results/psmc-test/psmc-plot-test/{sample}.eps",
        par="results/psmc-test/psmc-plot-test/{sample}.par"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc-test/psmc-plot-test/{sample}.log"
    benchmark:
        "results/benchmarks/psmc-test/psmc-plot-test/{sample}.bmk"
    shell:
        "psmc_plot.pl -u 8.0e-09 -g 3 {output.eps} {input} 2> {log}"


## rule to plot all PSMC outputs together!
# can also set min (-x) and max (-X) generations for mapping
rule psmc_plot_all_test:
    input:
        psmc=expand("results/psmc/run-psmc/{s}.psmc", s=sample_list),
    params:
        samps=expand("{s}, ", s=sample_list),
    output:
        "results/psmc-test/psmc-plot-all-test/all-together",
        #par="results/psmc-test/psmc-plot-all-test/all-together.par"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc-test/psmc-plot-test-test/all-together.log"
    benchmark:
        "results/benchmarks/psmc-test/psmc-plot-test-test/all-together.bmk"
    shell:
        " echo(psmc_plot.pl -u 8.0e-09 -g 3 -x 1000 -X 5000000 -M {params.samps} {input.psmc}) > {output} 2> {log}"