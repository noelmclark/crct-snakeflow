## rules to run admixture

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
rule fix_admixture_chroms:
    input:
        "results/plink/bed/aut-bisnps-no5indel.bim",
    output:
        flag="results/plink/bed/fix-chrom-flag.txt",
    params:
        pfx="results/plink/bed/aut-bisnps-no5indel"
    log:
        "results/logs/admixture/aut-bisnps-no5indel-fix-chrom.log"
    benchmark:
        "results/logs/admixture/aut-bisnps-no5indel-fix-chrom.bmk"
    shell:
        """
        ( mv {input} {input}.tmp && 
        awk '{{$1="0";print $0}}' {input}.tmp > {params.pfx}.bim && 
        rm {input}.tmp && 
        echo "admixture chroms fixed" > {output.flag} 
        ) 2> {log} 
        """


## runs through each of the selected k options -- super contrived
# admixture is weird and will not let you redirect the Q and P outputs - they will be produced in the current WD
# so we have to cd into where we want them to go. Also, the default for CV is 5 but can do up to 10
# input.flag makes sure the previous rule is run before trying this
# empty is so the directory gets made before we cd into it, and so the next rule knows to run this one first
rule test_k:
    input:
        bed="results/plink/bed/aut-bisnps-no5indel.bed",
        flag="results/plink/bed/fix-chrom-flag.txt",
    output:
        empty="results/admixture/CV_5/aut-bisnps-no5indel-{kclusters}.out",
    params:
        dir="results/admixture/CV_5/",
        pfx="aut-bisnps-no5indel-{kclusters}.out",
    conda:
        "../envs/admixture.yaml"
    resources:
        mem_mb=112200,
        time="23:59:59"
    threads:
        4
    log:
        "results/logs/admixture/CV_5/aut-snps-0.05-pruned-{kclusters}.log"
    benchmark:
        "results/benchmarks/admixture/CV_5/aut-snps-0.05-pruned-{kclusters}.bmk"
    shell:
        " ( > {output.empty} && "
        " cd {params.dir} && "
        " admixture --cv ../../../{input.bed} {wildcards.kclusters} -j{threads}> {params.pfx} ) 2> {log} " 

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


## Don't really understand how to use this at the moment but it works in an interative R session
## Rule to use Jonah Meier's admxiture plotting r script (https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r)
# currently works but not very pretty
# -i gives the info file that connects the sample names to pop/spp in the same order as admixture processed them
# -l specifies a comma-separated list of populations/species in the order you want them to be plotted
rule plot_admixture:
    input:
        pfx="results/admixture/CV_5/aut-bisnps-no5indel-",
        list="config/admixture-info.tsv"
    output:
        "results/admixture/rplot/aut-bisnps-no5indel"
    envmodules: 
        "R/4.2.2"
    log:
        "results/logs/admixture/rplot/aut-bisnps-no5indel.log"
    benchmark:
        "results/benchmarks/admixture/rplot/aut-bisnps-no5indel.bmk"
    shell:
        " Rscript workflow/scripts/admixture/plotADMIXTURE.r -p {input.pfx} "
        " -i {input.list} -k 3 -m 2 "
        " -l navajo,williamson,nanita,w_fk_boulder,dry_gulch,e_fk_piedra,steelman,abrams,roan,hunter,kelso,s_twin,w_antelope,severy,como,s_hayden,greenback,rio_grande,san_juan,bonneville,yellowstone,lahontan,coastal,westslope "
        " -o {output} 2> {log} "


### MAC1 ###
# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
rule fix_admixture_chroms_mac1:
    input:
        "results/plink/bed/MAC1/aut-bisnps-no5indel-nooutlier-MAC1.bim",
    output:
        flag="results/plink/bed/fix-chrom-flag-MAC1.txt",
    params:
        pfx="results/plink/bed/MAC1/aut-bisnps-no5indel-nooutlier-MAC1"
    log:
        "results/logs/admixture/aut-bisnps-no5indel-nooutlier-fix-chrom-MAC1.log"
    benchmark:
        "results/logs/admixture/aut-bisnps-no5indel-nooutlier-fix-chrom-MAC1.bmk"
    shell:
        """
        ( mv {input} {input}.tmp && 
        awk '{{$1="0";print $0}}' {input}.tmp > {params.pfx}.bim && 
        rm {input}.tmp && 
        echo "admixture chroms fixed" > {output.flag} 
        ) 2> {log} 
        """

rule test_k_mac1:
    input:
        bed="results/plink/bed/MAC1/aut-bisnps-no5indel-nooutlier-MAC1.bed",
        flag="results/plink/bed/fix-chrom-flag-MAC1.txt",
    output:
        empty="results/admixture/CV_5/MAC1/aut-bisnps-no5indel-nooutlier-MAC1-{kclusters}.out",
    params:
        dir="results/admixture/CV_5/MAC1/",
        pfx="aut-bisnps-no5indel-nooutlier-MAC1-{kclusters}.out",
    conda:
        "../envs/admixture.yaml"
    resources:
        mem_mb=112200,
        time="23:59:59"
    threads:
        4
    log:
        "results/logs/admixture/CV_5/MAC1/aut-bisnps-no5indel-nooutlier-MAC1-{kclusters}.log"
    benchmark:
        "results/benchmarks/admixture/CV_5/MAC1/aut-bisnps-no5indel-nooutlier-MAC1-{kclusters}.bmk"
    shell:
        " ( > {output.empty} && "
        " cd {params.dir} && "
        " admixture --cv ../../../../{input.bed} {wildcards.kclusters} -j{threads}> {params.pfx} ) 2> {log} " 

# grep-h CV log*.out to view cv values
rule get_best_k_mac1:
    input:
        expand("results/admixture/CV_5/MAC1/aut-bisnps-no5indel-nooutlier-MAC1-{k}", k=kclusters)
    output:
        "results/admixture/CV_5/MAC1/aut-bisnps-no5indel-nooutlier-MAC1.cv5.error"
    log:
        "results/logs/admixture/CV_5/MAC1/aut-bisnps-no5indel-nooutlier-MAC1.cv5.error.log"
    benchmark:
        "results/benchmarks/admixture/CV_5/MAC1/aut-bisnps-no5indel-nooutlier-MAC1.cv5.error.bmk"
    shell:
        " awk '/CV/ {print $3,$4}' {input} > {output} 2> {log} "

###############

### MAC3 ###
# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
rule fix_admixture_chroms_mac3:
    input:
        "results/plink/bed/MAC3/aut-bisnps-no5indel-nooutlier-MAC3.bim",
    output:
        flag="results/plink/bed/fix-chrom-flag-MAC3.txt",
    params:
        pfx="results/plink/bed/MAC3/aut-bisnps-no5indel-nooutlier-MAC3"
    log:
        "results/logs/admixture/aut-bisnps-no5indel-nooutlier-fix-chrom-MAC3.log"
    benchmark:
        "results/logs/admixture/aut-bisnps-no5indel-nooutlier-fix-chrom-MAC3.bmk"
    shell:
        """
        ( mv {input} {input}.tmp && 
        awk '{{$1="0";print $0}}' {input}.tmp > {params.pfx}.bim && 
        rm {input}.tmp && 
        echo "admixture chroms fixed" > {output.flag} 
        ) 2> {log} 
        """

rule test_k_mac3:
    input:
        bed="results/plink/bed/MAC3/aut-bisnps-no5indel-nooutlier-MAC3.bed",
        flag="results/plink/bed/fix-chrom-flag-MAC3.txt",
    output:
        empty="results/admixture/CV_5/MAC3/aut-bisnps-no5indel-nooutlier-MAC3-{kclusters}.out",
    params:
        dir="results/admixture/CV_5/MAC3/",
        pfx="aut-bisnps-no5indel-nooutlier-MAC3-{kclusters}.out",
    conda:
        "../envs/admixture.yaml"
    resources:
        mem_mb=112200,
        time="23:59:59"
    threads:
        4
    log:
        "results/logs/admixture/CV_5/MAC3/aut-bisnps-no5indel-nooutlier-MAC3-{kclusters}.log"
    benchmark:
        "results/benchmarks/admixture/CV_5/MAC3/aut-bisnps-no5indel-nooutlier-MAC3-{kclusters}.bmk"
    shell:
        " ( > {output.empty} && "
        " cd {params.dir} && "
        " admixture --cv ../../../../{input.bed} {wildcards.kclusters} -j{threads}> {params.pfx} ) 2> {log} " 

# grep-h CV log*.out to view cv values
rule get_best_k_mac3:
    input:
        expand("results/admixture/CV_5/MAC3/aut-bisnps-no5indel-nooutlier-MAC3-{k}", k=kclusters)
    output:
        "results/admixture/CV_5/MAC3/aut-bisnps-no5indel-nooutlier-MAC3.cv5.error"
    log:
        "results/logs/admixture/CV_5/MAC3/aut-bisnps-no5indel-nooutlier-MAC3.cv5.error.log"
    benchmark:
        "results/benchmarks/admixture/CV_5/MAC3/aut-bisnps-no5indel-nooutlier-MAC3.cv5.error.bmk"
    shell:
        " awk '/CV/ {print $3,$4}' {input} > {output} 2> {log} "

###############

### MAC5 ###
# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
rule fix_admixture_chroms_mac5:
    input:
        "results/plink/bed/MAC5/aut-bisnps-no5indel-nooutlier-MAC5.bim",
    output:
        flag="results/plink/bed/fix-chrom-flag-MAC5.txt",
    params:
        pfx="results/plink/bed/MAC5/aut-bisnps-no5indel-nooutlier-MAC5"
    log:
        "results/logs/admixture/aut-bisnps-no5indel-nooutlier-fix-chrom-MAC5.log"
    benchmark:
        "results/logs/admixture/aut-bisnps-no5indel-nooutlier-fix-chrom-MAC5.bmk"
    shell:
        """
        ( mv {input} {input}.tmp && 
        awk '{{$1="0";print $0}}' {input}.tmp > {params.pfx}.bim && 
        rm {input}.tmp && 
        echo "admixture chroms fixed" > {output.flag} 
        ) 2> {log} 
        """

rule test_k_mac5:
    input:
        bed="results/plink/bed/MAC5/aut-bisnps-no5indel-nooutlier-MAC5.bed",
        flag="results/plink/bed/fix-chrom-flag-MAC5.txt",
    output:
        empty="results/admixture/CV_5/MAC5/aut-bisnps-no5indel-nooutlier-MAC5-{kclusters}.out",
    params:
        dir="results/admixture/CV_5/MAC5/",
        pfx="aut-bisnps-no5indel-nooutlier-MAC5-{kclusters}.out",
    conda:
        "../envs/admixture.yaml"
    resources:
        mem_mb=112200,
        time="23:59:59"
    threads:
        4
    log:
        "results/logs/admixture/CV_5/MAC5/aut-bisnps-no5indel-nooutlier-MAC5-{kclusters}.log"
    benchmark:
        "results/benchmarks/admixture/CV_5/MAC5/aut-bisnps-no5indel-nooutlier-MAC5-{kclusters}.bmk"
    shell:
        " ( > {output.empty} && "
        " cd {params.dir} && "
        " admixture --cv ../../../../{input.bed} {wildcards.kclusters} -j{threads}> {params.pfx} ) 2> {log} " 

# grep-h CV log*.out to view cv values
#" awk '/CV/ {print $3,$4}' {input} > {output} 2> {log} "
rule get_best_k_mac5:
    input:
        expand("results/admixture/CV_5/MAC5/aut-bisnps-no5indel-nooutlier-MAC5-{k}.out", k=kclusters)
    output:
        "results/admixture/CV_5/aut-bisnps-no5indel-nooutlier-MAC5.cv5.error"
    log:
        "results/logs/admixture/CV_5/aut-bisnps-no5indel-nooutlier-MAC5.cv5.error.log"
    benchmark:
        "results/benchmarks/admixture/CV_5/aut-bisnps-no5indel-nooutlier-MAC5.cv5.error.bmk"
    shell:
        " grep-h CV {input} > output 2> {log} "

###############