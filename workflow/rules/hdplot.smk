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
        bcftools query -f '[%AD\\t]\\n' {input.bcf} > {output.tsv} 2> {log}
        """

rule compute_hdplot_stats:
    input:
        genotypes="results/hdplot/all-samples/all-genotypes.tsv",
        depths="results/hdplot/all-samples/all-allele-depths.tsv"
    output:
        table="results/hdplot/all-samples/all-hdplot-output.tsv"
    log:
        "results/logs/hdplot/process-hdplot/all-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/hdplot/process-hdplot/all-aut-bisnps-no5indel.bmk",
    script:
        " workflow/scripts/hdplot/process_hdplot.awk {input.genotypes} {input.depths} > {output.table} 2> {log} "

###############################################################################################