## this script implements a version of HDplot using bcftools commands
## based on the original at https://github.com/gjmckinney/HDplot/tree/master
## and written in this paper: https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12613 

# the top set runs on the all sample bcf

rule extract_all_genotypes:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
    output:
        tsv="results/hdplot/all-samples/all-genotypes.tsv"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hdplot/extract-genotypes/all-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/hdplot/extract-genotypes/all-aut-bisnps-no5indel.bmk",
    shell:
        """
        bcftools query -f '%CHROM\\t%POS\\t%ID[\\t%GT]\\n' {input.bcf} > {output.tsv} 2> {log}
        """

rule extract_all_allele_depths:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
    output:
        tsv="results/hdplot/all-samples/all-allele-depths.tsv"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hdplot/extract-allele-depths/all-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/hdplot/extract-allele-depths/all-aut-bisnps-no5indel.bmk",
    shell:
        """
        bcftools query -f '[%AD\\t]\\n' {input.bcf} | sed 's/\\t$//' > {output.tsv} 2> {log}
        """

rule compute_parallel_hdplot_stats:
    input:
        geno="results/hdplot/all-samples/all-genotypes.tsv",
        depths="results/hdplot/all-samples/all-allele-depths.tsv"
    output:
        "results/hdplot/all-samples/all-hdplot-output.tsv"
    resources:
        mem_mb=90000,
    log:
        "results/logs/hdplot/process-hdplot/all-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/hdplot/process-hdplot/all-aut-bisnps-no5indel.bmk",
    shell:
        " awk -f workflow/scripts/hdplot/hdplot_parallel.awk {input.geno} {input.depths} > {output} 2> {log} "

###############################################################################################

# the bottom set runs on the bcf filtered to only yellowstone group alleles -- needs editing

rule filter_hdplot_srm_bcf:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
        samps="config/hdplot-srm-samps.txt",
    output:
        bcf="results/bcf/hdplot-srm-aut-bisnp-no5indel.bcf",
        tbi="results/bcf/hdplot-srm-aut-bisnp-no5indel.bcf.csi",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hdplot/filter-srm-bcf/srm-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/hdplot/filter-srm-bcf/srm-aut-bisnps-no5indel.bmk",
    shell:
        """
        bcftools view -Ou -S {input.samps} {input.bcf} |
        bcftools +fill-tags -Ob -- -t all > {output.bcf} && 
        bcftools index {output.bcf} 2> {log}
        """

rule extract_srm_genotypes:
    input:
        bcf="results/bcf/hdplot-srm-aut-bisnp-no5indel.bcf",
        tbi="results/bcf/hdplot-srm-aut-bisnp-no5indel.bcf.csi",
    output:
        tsv="results/hdplot/srm-samples/srm-genotypes.tsv"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hdplot/extract-genotypes/srm-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/hdplot/extract-genotypes/srm-aut-bisnps-no5indel.bmk",
    shell:
        """
        bcftools query -f '%CHROM\\t%POS\\t%ID[\\t%GT]\\n' {input.bcf} > {output.tsv} 2> {log}
        """

rule extract_srm_allele_depths:
    input:
        bcf="results/bcf/hdplot-srm-aut-bisnp-no5indel.bcf",
        tbi="results/bcf/hdplot-srm-aut-bisnp-no5indel.bcf.csi",
    output:
        tsv="results/hdplot/srm-samples/srm-allele-depths.tsv"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hdplot/extract-allele-depths/srm-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/hdplot/extract-allele-depths/srm-aut-bisnps-no5indel.bmk",
    shell:
        """
        bcftools query -f '[%AD\\t]\\n' {input.bcf} | sed 's/\\t$//' > {output.tsv} 2> {log}
        """

rule compute_srm_parallel_hdplot_stats:
    input:
        geno="results/hdplot/srm-samples/srm-genotypes.tsv",
        depths="results/hdplot/srm-samples/srm-allele-depths.tsv"
    output:
        "results/hdplot/srm-samples/srm-hdplot-output.tsv"
    resources:
        mem_mb=90000,
    log:
        "results/logs/hdplot/process-hdplot/srm-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/hdplot/process-hdplot/srm-aut-bisnps-no5indel.bmk",
    shell:
        " awk -f workflow/scripts/hdplot/hdplot_parallel.awk {input.geno} {input.depths} > {output} 2> {log} "

###############################################################################################