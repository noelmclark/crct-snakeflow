## following rules are to run PSMC
# based on lh3 documentation at: https://github.com/lh3/psmc

## rule to get a consensus fastq sequence file for PSMC
# option -C 50 downgrades mapping quality (by coeff given) for reads containing excessive mismatches
# option -d sets and minimum read depth and -D sets the maximum 
# It is recommended to set -d to a third of the average depth and -D to twice
# this takes the bams from each individual and generates a tmp vcf file then 
# generates a consensus sequence for psmc
# alternatively add a rule that does this from the hard filtered vcf files after
# joint calling them in the calling phase
rule psmc_consensus_sequence:
    input:
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",  
        ref="resources/genome/OmykA.fasta",
    output:
        temp("results/psmc/psmc-consensus-sequence/{sample}.fq.gz") #temp tells snakemake to remove these files once they're no longer needed downstream
    conda:
        "../envs/sambcftools.yaml"
    resources:
        time="23:59:59"
    log:
        "results/logs/psmc/psmc-consensus-sequence/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/psmc-consensus-sequence/{sample}.bmk"
    shell:
        "samtools mpileup -C50 -uf {input.ref} {input.bam} | bcftools call -c - | " 
        "vcfutils.pl vcf2fq -d 6 -D 36 | gzip > {output} 2> {log}"


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


# rule to run psmc
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
        "psmc -N25 -t15 -r5 -p '4+25*2+4+6' -o {output} {input} 2> {log}"


## rule to plot psmc to visualize result
# -u [per-generation mutation rate] from https://doi.org/10.1371/journal.pgen.1010918
# -g [generation time in years] 
rule psmc_plot:
    input:
        "results/psmc/run-psmc/{sample}.psmc"
    output:
        eps="results/psmc/psmc-plot/{sample}.eps",
        par="results/psmc/psmc-plot/{sample}.par"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/psmc-plot/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/psmc-plot/{sample}.bmk"
    shell:
        "psmc_plot.pl -u 8.0e-09 -g 3 {output.eps} {input} 2> {log}"


## alternative rule to run psmc2history & history2ms to generate the ms
# command line that simulates the history inferred by PSMC
# I use the rule psmc_plot instead
rule psmc2history2ms:
    input:
        "results/psmc/run-psmc/{sample}.psmc"
    output:
        "results/psmc/psmc2history2ms/{sample}-ms-cmd.sh"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/psmc2history2ms/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/psmc2history2ms/{sample}.bmk"
    shell:
        "psmc2history.pl {input} | history2ms.pl > {output} 2> {log}"