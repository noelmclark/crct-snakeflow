

# draft rule for vcf 2 psmc consensus sequence vs from the bams 
#rule vcf2_psmc_consensus_sequence:
#    input:
#        vcf="results/
#    shell:
#        "{input.ref} {input.bam}" 
#        "vcfutils.pl vcf2fq -d 6 -D 36 | gzip > {output} 2> {log}"  

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





# rule to run psmc2history & history2ms to generate the ms
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





# rule to plot psmc to visualize result
# -u [per-generation mutation rate] from https://doi.org/10.1371/journal.pgen.1010918
# -g [generation time in years] 
rule psmc_plot:
    input:
        "results/psmc/run-psmc/{sample}.psmc"
    output:
        "results/psmc/psmc-plot/{sample}.eps"
    conda:
        "../envs/psmc.yaml"
    log:
        "results/logs/psmc/psmc-plot/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/psmc-plot/{sample}.bmk"
    shell:
        "psmc_plot.pl -u 8.0e-09 -g 3 {output} {input} 2> {log}"