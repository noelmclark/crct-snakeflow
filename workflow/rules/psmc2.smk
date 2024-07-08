

# draft rule for vcf 2 psmc consensus sequence vs from the bams 
#rule vcf2_psmc_consensus_sequence:
#    input:
#        vcf="results/
#    shell:
#        "{input.ref} {input.bam}" 
#        "vcfutils.pl vcf2fq -d 6 -D 36 | gzip > {output} 2> {log}"  



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