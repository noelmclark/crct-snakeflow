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
    envmodules: "R/4.4.0"
    output:
        tsv="results/rCNV/detect_deviants_out.tsv"
    log:
    	"results/logs/rCNV/detect_deviants.log"
    script:
    	"../scripts/rCNV.R"