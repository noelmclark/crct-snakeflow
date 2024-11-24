## rules to run admixture

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
rule fix_admixture_chroms:
    input:
        "results/plink/bed/aut-bisnps-no5indel.bim"
    output:
        pfx="results/plink/bed/aut-bisnps-no5indel",
        flag="results/plink/bed/fix-chrom-flag.txt"
    log:
        "results/logs/admixture/aut-bisnps-no5indel-fix-chrom.log"
    benchmark:
        "results/logs/admixture/aut-bisnps-no5indel-fix-chrom.bmk"
    shell:
        " ( mv {input} {input}.tmp && "
        " awk '{{$1=/"0/";print $0}}' {input}.tmp > {output.pfx}.bim && "
        " rm {input}.tmp && "
        " echo /"admixture chroms fixed/" > {output.flag}) 2> {log} "


## runs through each of the selected k options
# admixture is weird and will not let you redirect the Q and P outputs - they will be produced in the current WD
# so we have to cd into where we want them to go. Also, the default for CV is 5 but can do up to 10
# I might have to reroute to the input file in the shell code after cd-ing into the output dir, we'll see
# input.flag makes sure the previous rule is run before trying this
# empty is so I can ask for the next rule and it knows to run this one first
rule test_k:
    input:
        bed="results/plink/bed/aut-bisnps-no5indel.bed",
        flag="results/plink/bed/fix-chrom-flag.txt"
    output:
        pfx="aut-bisnps-no5indel-{kclusters}.out",
        empty="results/admixture/CV_5/aut-bisnps-no5indel-{kclusters}.out"
    params:
        dir="results/admixture/CV_5/",
    conda:
        "../envs/admixture.yaml"
    resources:
        mem_mb=11220,
        cpus=2,
        time="12:00:00"
    log:
        "results/logs/admixture/CV_5/aut-snps-0.05-pruned-{kclusters}.log"
    benchmark:
        "results/benchmarks/admixture/CV_5/aut-snps-0.05-pruned-{kclusters}.bmk"
    shell:
        " cd {params.dir} && "
        " admixture --cv {input.bed} {wildcards.kclusters} > {params.dir}{output.pfx} 2> {log} " 

# grep-h CV log*.out to view cv values
rule get_best_k:
    input:
        expand("results/admixture/CV_5/aut-snps-0.05-pruned{k}.out", k=kclusters)
    output:
        "results/admixture/CV_5/aut-snps-0.05-pruned.cv5.error"
    log:
        "results/logs/admixture/CV_5/aut-snps-0.05-pruned.cv5.error.log"
    benchmark:
        "results/benchmarks/admixture/CV_5/aut-snps-0.05-pruned.cv5.error.bmk"
    shell:
        " awk '/CV/ {print $3,$4}' {input} > {output} 2> {log} "

## Rule to use Jonah Meier's admxiture plotting r script (https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r)
# -l specifies a comma-separated list of populations/species in the order to be plotted
#rule plot_admixture:
#    input:
#        pfx="results/admixture/aut-snps-0.05-pruned-"
#        list="admixture-info.tsv"
#    output:
#        "results/admixture/rplot/aut-snps-0.05-pruned"
#    envmodules: 
#        "R/4.2.2"
#    log:
#    benchmark:
#    shell:
#        " /scripts/admixture/plotADMIXTURE.r -p {input.pfx} -i {input.list} -k 24 -m 6 -l ?? -o {output}"