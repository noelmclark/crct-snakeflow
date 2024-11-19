## rules to run admixture

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
#rule fix_admixture_chroms:
#    input:
#        "results/plink/pca/aut-snps-0.05-pruned-pca.bim"
#    output:
#        "results/plink/pca/aut-snps-0.05-pruned-pca"
#    log:
#        "results/logs/admixture/aut-snps-0.05-pruned-pca-fix-chrom.log"
#    benchmark:
#        "results/logs/admixture/aut-snps-0.05-pruned-pca-fix-chrom.bmk"
#    shell:
#        " ( mv {input} {input}.tmp && "
#        " awk '{$1="0";print $0}' {input}.tmp > {output}.bim && "
#        " rm {input}.tmp ) 2> {log} "


#rule admixture_first_k:
#    input:
#        "results/plink/pca/aut-snps-0.05-pruned-pca.bed"
#    output:
#        "results/admixture/aut-snps-0.05-pruned-pca.10.Q",
#    conda:
#        "../envs/admixture.yaml"
#    log:
#        "results/logs/admixture/aut-snps-0.05-pruned-pca.10.log"
#    benchmark:
#        "results/benchmarks/admixture/aut-snps-0.05-pruned-pca.10.bmk"
#    shell:
#        "admixture --cv {input} 10"

## runs through each of the selected k options
# grep-h CV log*.out to view cv values
rule test_k:
    input:
        "results/plink/bed/aut-snps-0.05-pruned.bed"
    output:
        "results/admixture/test_k/aut-snps-0.05-pruned"
    conda:
        "../envs/admixture.yaml"
    log:
        "results/logs/admixture/test_k/aut-snps-0.05-pruned.log"
    benchmark:
        "results/benchmarks/admixture/test_k/aut-snps-0.05-pruned.bmk"
    shell:
        " (for k in {wildcards.kclusters}; do "
        " admixture --cv {input} $k > {output}$k.out "
        " done ) 2> {log} " 

#rule get_best_k:
    input:
        expand("results/admixture/test_k/aut-snps-0.05-pruned{k}.out", k=kclusters)
    output:
        "results/admixture/test_k/aut-snps-0.05-pruned.cv5.error"
    log:
        "results/logs/admixture/test_k/aut-snps-0.05-pruned.cv5.error.log"
    benchmark:
        "results/benchmarks/admixture/test_k/aut-snps-0.05-pruned.cv5.error.bmk"
    shell:
        " awk '/CV/ {print $3,$4}' {input}* | cut -c 4,7-20 > {output} 2> {log} "

#rule plot_admixture:
#    input:
#    output:
#    envmodules: 
#        "R/4.2.2"
#    log:
#    benchmark:
#    shell:
#        " /scripts/admixture/plotADMIXTURE.r -p $FILE -i $FILE.list -k 5 -l "