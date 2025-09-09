## going to try running rCNV on the BCF split by scatters (scatter_idx) which is
## based on Eric's post-bcf workflow https://github.com/eriqande/mega-post-bcf-exploratory-snakeflows
## first rule creates the regions files in the format that bcftools needs to split up the BCF 
## second rule splits the given BCF into scatter index regions of 120000bp length
## third rule is not used currently because I run rCNV in an acompile node using the rCNV.R script and R in a conda env

rule make_rCNV_scat_regions:
    input:
        scat_path="results/scatter_config/scatters_100000.tsv"
    params:
        scat="{scatter}"
    output:
        "results/rCNV-by-scat/scat-regions/{scatter}.scat-regions.tsv",
    log:
        "results/logs/rCNV-by-scat/scat-regions/{scatter}.scat-regions.log"
    benchmark:
        "results/benchmarks/rCNV-by-scat/scat-regions/{scatter}.scat-regions.bmk"
    shell:
        " (awk -v scat='{params.scat}' -f workflow/scripts/pca/get_scat_regions.awk {input.scat_path} > {output}) 2> {log} "


rule rCNV_get_scatters:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
        regions="results/rCNV-by-scat/scat-regions/{scatter}.scat-regions.tsv",
    output:
        vcf="results/rCNV-by-scat/vcf/aut-bisnp-no5indel-{scatter}.vcf",
    log:
        "results/logs/rCNV-by-scat/vcf/aut-bisnp-no5indel-{scatter}.log",
    benchmark:
        "results/benchmarks/rCNV-by-scat/vcf/aut-bisnp-no5indel-{scatter}.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools view -Ov -R {input.regions} {input.bcf} > {output.vcf} 2> {log} "


## HERE is where I go into an interactive remote R session on Alpine to run ##
## the rCNV.R script on each of the rCNV scatter VCFs from the previous rule ##
## to generate the allele info file needed for the next rule ##

rule rCNV_dvs_cnv:
    input:
        vcf="results/rCNV-by-scat/vcf/aut-bisnp-no5indel-{scatter}.vcf",
        ainfo="results/rCNV-by-scat/out/aut-bisnp-no5indel-{scatter}_allele_info.tsv",
    envmodules: 
        "R/4.2.2"
    output:
        dvs="results/rCNV-by-scat/dvs/aut-bisnp-no5indel-{scatter}_dvs.tsv",
        cnv="results/rCNV-by-scat/cnv/aut-bisnp-no5indel-{scatter}_cnv.tsv",
    log:
    	"results/logs/rCNV-by-scat/dvs-cnv/dvs-cnv-{scatter}.log",
    benchmark:
        "results/benchmarks/rCNV-by-scat/dvs-cnv/dvs-cnv-{scatter}.bmk",
    script:
    	"../scripts/rCNV/rCNV-dvs-cnv.R"