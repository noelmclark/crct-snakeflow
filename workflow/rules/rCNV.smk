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
## to generate the allele info file needed for the next script ##
## then I use the rCNV-dvs-cnv.R file to ID deviants/cnvs and clean them ##

rule rCNV_dvs_cnv_targets:
    input:
        dvs="results/rCNV-by-scat/dvs-clean-merged/aut-bisnp-no5indel_dvs_clean_merged.tsv",
        cnv="results/rCNV-by-scat/cnv-clean-merged/aut-bisnp-no5indel_cnv_clean_merged.tsv",
    output:
        dvs="results/rCNV-by-scat/dvs-cnv-targets/aut-bisnp-no5indel_dvs_targets.tsv",
        cnv="results/rCNV-by-scat/dvs-cnv-targets/aut-bisnp-no5indel_cnv_targets.tsv",
    log:
        "results/logs/rCNV-by-scat/dvs-cnv-targets/dvs-cnv-targets.log",
    benchmark:
        "results/benchmarks/rCNV-by-scat/dvs-cnv-targets/dvs-cnv-targets.bmk",
    shell:
        "((awk -f workflow/scripts/rCNV/get-rCNV-targets.awk {input.cnv} | sort -k1,1 -k2,2n > {output.cnv}); "
        "(awk -f workflow/scripts/rCNV/get-rCNV-targets.awk {input.dvs} | sort -k1,1 -k2,2n > {output.dvs})) 2> {log} "


rule rCNV_index_targets:
    input:
        dvs="results/rCNV-by-scat/dvs-cnv-targets/aut-bisnp-no5indel_dvs_targets.tsv",
        cnv="results/rCNV-by-scat/dvs-cnv-targets/aut-bisnp-no5indel_cnv_targets.tsv",
    output:
        dvs="results/rCNV-by-scat/dvs-cnv-targets/aut-bisnp-no5indel_dvs_targets.tsv.gz",
        cnv="results/rCNV-by-scat/dvs-cnv-targets/aut-bisnp-no5indel_cnv_targets.tsv.gz",
    conda:
        "../envs/sambcftools.yaml",
    log:
        "results/logs/rCNV-by-scat/dvs-cnv-targets/dvs-cnv-index-targets.log",
    benchmark:
        "results/benchmarks/rCNV-by-scat/dvs-cnv-targets/dvs-cnv-index-targets.bmk",
    shell:
        "((bgzip -c {input.dvs} > {output.dvs} && tabix -s1 -b2 -e2 {output.dvs} ); "
        "(bgzip -c {input.cnv} > {output.cnv} && tabix -s1 -b2 -e2 {output.cnv} )) 2> {log} "

rule rCNV_filter_bcf_by_dvs:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
        targets="results/rCNV-by-scat/dvs-cnv-targets/aut-bisnp-no5indel_dvs_targets.tsv.gz",
    output:
        bcf="results/bcf/aut-bisnps-no5indel-rcnv-by-dvs.bcf",
    log:
        "results/logs/rCNV-by-scat/bcf/aut-bisnp-no5indel-rcnv-by-dvs.log",
    benchmark:
        "results/benchmarks/rCNV-by-scat/bcf/aut-bisnp-no5indel-rcnv-by-dvs.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools filter -Ob -T {input.targets} {input.bcf} > {output.bcf}:
        " bcftools index {output.bcf}) 2> {log} "

rule rCNV_filter_bcf_by_cnv:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
        targets="results/rCNV-by-scat/dvs-cnv-targets/aut-bisnp-no5indel_cnv_targets.tsv.gz",
    output:
        bcf="results/bcf/aut-bisnps-no5indel-rcnv-by-cnv.bcf",
    log:
        "results/logs/rCNV-by-scat/bcf/aut-bisnp-no5indel-rcnv-by-cnv.log",
    benchmark:
        "results/benchmarks/rCNV-by-scat/bcf/aut-bisnp-no5indel-rcnv-by-cnv.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools filter -Ob -T {input.targets} {input.bcf} > {output.bcf}:
        " bcftools index {output.bcf}) 2> {log} "