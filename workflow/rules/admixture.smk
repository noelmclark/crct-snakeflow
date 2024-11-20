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
# admixture is weird and will not let you redirect the Q and P outputs - they will be produced in the current WD
# so we have to cd into where we want them to go. Also, the default for CV is 5 but can do up to 10

# right now I just run these locally in a compute node
rule test_k:
    input:
        bed="results/plink/bed/aut-snps-0.05-pruned.bed",
        empty="results/admixture/test_k/CV_5/aut-snps-0.05-pruned-{kclusters}.out"
    output:
        pfx="aut-snps-0.05-pruned-{kclusters}"
    params:
        dir="results/admixture/test_k/CV_5/",
    conda:
        "../envs/admixture.yaml"
    resources:
        mem_mb=9400,
        cpus=2,
    log:
        "results/logs/admixture/test_k/CV_5/aut-snps-0.05-pruned-{kclusters}.log"
    benchmark:
        "results/benchmarks/admixture/test_k/CV_5/aut-snps-0.05-pruned-{kclusters}.bmk"
    shell:
        " cd {params.dir} && "
        " admixture --cv {input.bed} {wildcards.kclusters} > {params.dir}{output.pfx}.out 2> {log} " 

# grep-h CV log*.out to view cv values
rule get_best_k:
    input:
        expand("results/admixture/test_k/CV_5/aut-snps-0.05-pruned{k}.out", k=kclusters)
    output:
        "results/admixture/test_k/CV_5/aut-snps-0.05-pruned.cv5.error"
    log:
        "results/logs/admixture/test_k/CV_5/aut-snps-0.05-pruned.cv5.error.log"
    benchmark:
        "results/benchmarks/admixture/test_k/CV_5/aut-snps-0.05-pruned.cv5.error.bmk"
    shell:
        " awk '/CV/ {print $3,$4}' {input} > {output} 2> {log} "

## Rule to use Jonah Meier's admxiture plotting r script (https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r)
# -l specifies a comma-separated list of populations/species in the order to be plotted
#rule plot_admixture:
#    input:
#        pfx="results/admixture/test_k/aut-snps-0.05-pruned-"
#        list="admixture-info.tsv"
#    output:
#        "results/admixture/rplot/aut-snps-0.05-pruned"
#    envmodules: 
#        "R/4.2.2"
#    log:
#    benchmark:
#    shell:
#        " /scripts/admixture/plotADMIXTURE.r -p {input.pfx} -i {input.list} -k 24 -m 6 -l ?? -o {output}"