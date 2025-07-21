# This uses rCNV described here: https://piyalkarum.github.io/rCNV/articles/rCNV.html
# Needs to have R installed and on the path.

rule rCNV_get_vcf:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
    output:
        vcf="results/bcf/aut-bisnps-no5indel.vcf",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/rCNV/get-vcf/all-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/rCNV/get-vcf/all-aut-bisnps-no5indel.bmk",
    shell:
        """
        bcftools view -Ov {input.bcf} > {output.vcf} 2> {log}
        """

rule rCNV_detect_deviants:
    input:
        vcf="results/bcf/aut-bisnps-no5indel.vcf",
    conda:
        "../envs/rCNV.yaml",
    envmodules: 
        "R/4.2.2"
    output:
        tsv="results/rCNV/detect-out.tsv"
    log:
    	"results/logs/rCNV/detect_deviants.log"
    script:
    	"../scripts/rCNV.R"

rule rCNV_test:
    input:
        vcf="results/bcf/aut-bisnps-no5indel.vcf",
    envmodules: 
        "R/4.2.2"
    output:
        tsv="results/rCNV/test_out.tsv"
    resources:
        mem_mb=90000,
    log:
    	"results/logs/rCNV/test.log"
    script:
    	"../scripts/rCNV-test.R"


##################################################################3

## going to try running rCNV on the BCF split by chromosome which I used in the hpsmc-chroms config file 
## based on Eric's post-bcf workflow https://github.com/eriqande/mega-post-bcf-exploratory-snakeflows
## first rule splits the given BCF into autosomal chroms
## second rule runs rCNV on each of these mini vcf files

rule rCNV_get_scatters:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
    output:
        vcf="results/rCNV-by-scat/vcf/aut-bisnp-no5indel-{hpsmcchroms}.vcf",
    log:
        "results/logs/rCNV-by-scat/vcf/aut-bisnp-no5indel-{hpsmcchroms}.log",
    benchmark:
        "results/benchmarks/rCNV-by-scat/vcf/aut-bisnp-no5indel-{hpsmcchroms}.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools view -Ov -r {wildcards.hpsmcchroms} {input.bcf} > {output.vcf} 2> {log} "


rule rCNV_test_run_by_scatters:
    input:
        vcf="results/rCNV-by-scat/vcf/aut-bisnp-no5indel-{hpsmcchroms}.vcf",
    envmodules: 
        "R/4.2.2"
    output:
        tsv="results/rCNV-by-scat/test-{hpsmcchroms}-deviants-out.tsv",
    log:
    	"results/logs/rCNV-by-scat/test-{hpsmcchroms}-deviants-out.log",
    benchmark:
        "results/benchmarks/rCNV-by-scat/test-{hpsmcchroms}-deviants-out.bmk",
    script:
    	"../scripts/rCNV-test.R"