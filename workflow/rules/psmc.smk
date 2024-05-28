# rule to get a consensus fastq sequence file for PSMC
# option -d sets and minimum read depth and -D sets the maximum 
# It is recommended to set -d to a third of the average depth and -D to twice
rule psmc_consensus_sequence:
    input:
        bam="results/mkdup/{sample}.bam",
        ref="resources/genome/OmykA.fasta",
    output:
        "results/psmc-consensus-sequence/{sample}.fq.gz"
    conda:
        "../envs/sambcftools.yaml"
    resources:
        time="23:59:59"
    log:
        "results/psmc-consensus-sequence/{sample}.log"
    benchmark:
        "results/benchmarks/psmc-consensus-sequence/{sample}.bmk"
    shell:
        "bcftools mpileup -C50 -f {input.ref} {input.bam} | bcftools call -c - | " 
        "vcfutils.pl vcf2fq -d 6 -D 36 | gzip > {output} 2> {log}"


# following rules are to run PSMC
# based on lh3 documentation at: https://github.com/lh3/psmc

# rule to create psmcfa file per sample
rule psmcfa:
    input:
        "results/psmc-consensus-sequence/{sample}.fq.gz"
    output:
        "results/psmcfa/{sample}.psmcfa"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/psmcfa/{sample}.log"
    benchmark:
        "results/benchmarks/psmcfa/{sample}.bmk"
    shell:
        "fq2psmcfa -q20 {input} > {output} 2> {log}"



# rule to run psmc
rule run_psmc:
    input:
        "results/psmcfa/{sample}.psmcfa"
    output:
        "results/run-psmc/{sample}.psmc"
    conda:
        "envs/psmc.yaml"
    log:
        "results/run-psmc/{sample}.log"
    benchmark:
        "results/benchmarks/run-psmc/{sample}.bmk"
    shell:
        "psmc -N25 -t15 -r5 -p '4+25*2+4+6' -o {output} {input} 2> {log}"





# rule to run psmc2history & history2ms to generate the ms
# command line that simulates the history inferred by PSMC
rule psmc2history2ms:
    input:
        "results/run-psmc/{sample}.psmc"
    output:
        "results/psmc2history2ms/{sample}-ms-cmd.sh"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/psmc2history2ms/{sample}.log"
    benchmark:
        "results/benchmarks/psmc2history2ms/{sample}.bmk"
  shell:
    "psmc2history.pl {input} | history2ms.pl > {output} 2> {log}"





# rule to plot psmc to visualize result
# -u [per-generation mutation rate] from https://doi.org/10.1371/journal.pgen.1010918
# -g [generation time in years] 
rule psmc_plot:
    input:
        "results/run-psmc/{sample}.psmc"
    output:
        "results/psmc-plot/{sample}.eps"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/psmc-plot/{sample}.log"
    benchmark:
        "results/benchmarks/psmc-plot/{sample}.bmk"
    shell:
        "psmc_plot.pl -u 8.0e-09 -g 3 {output} {input} 2> {log}"