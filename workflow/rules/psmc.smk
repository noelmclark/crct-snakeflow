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
        bam="results/angsd_bams/overlap_clipped/{sample}.bam", #should this be a bam with -baq 2 done on it for indel stuff? 
        ref="resources/genome/OmykA.fasta",
    output:
        "results/psmc/bams2psmc/psmc-consensus-sequence/{sample}.fq.gz"
    conda:
        "../envs/sambcftools.yaml"
    resources:
        time="23:59:59"
    log:
        "results/logs/psmc/bams2psmc/psmc-consensus-sequence/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/bams2psmc/psmc-consensus-sequence/{sample}.bmk"
    shell:
        "bcftools mpileup -C50 -f {input.ref} {input.bam} | bcftools call -c - | " 
        "vcfutils.pl vcf2fq -d 6 -D 36 | gzip > {output} 2> {log}"


# rule to create psmcfa file per sample
rule psmcfa:
    input:
        "results/psmc/bams2psmc/psmc-consensus-sequence/{sample}.fq.gz"
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
        "envs/psmc.yaml"
    log:
        "results/logs/psmc/run-psmc/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/run-psmc/{sample}.bmk"
    shell:
        "psmc -N25 -t15 -r5 -p '4+25*2+4+6' -o {output} {input} 2> {log}"


## rule to run psmc2history & history2ms to generate the ms
# command line that simulates the history inferred by PSMC
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