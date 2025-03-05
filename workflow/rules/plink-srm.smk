#### This rule file is used to subset our original BCF (aut-bisnp-no5indel.bcf) to just
### samples from the SRM and their associated variants - identified as necessary following defense

# this rule generates a global allele freq file
# that is necessary for running PCA with all samples 
rule make_srm_pgen:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
        srmfile="config/plink-srm-samples.tsv",
    output:
        pgen="results/plink/srm-subset/make-pgen/srm-aut-bisnps-no5indel",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=7480    
    log:
        "results/logs/plink/srm-subset/make-pgen/srm-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/plink/srm-subset/make-pgen/srm-aut-bisnps-no5indel.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --keep {input.srmfile} "
        " --make-pgen "
        " --out {output.pgen} 2> {log} "

rule srm_filter_missingness:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
        srmfile="config/plink-srm-samples.tsv",
    output:
        acount="results/plink/srm-subset/allele-count/srm-aut-bisnps-no5indel",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=7480    
    log:
        "results/logs/plink/srm-subset/allele-count/srm-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/plink/srm-subset/allele-count/srm-aut-bisnps-no5indel.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --keep {input.srmfile} "
        " --geno 0.1 "
        " --freq counts "
        " --out {output.acount} 2> {log} "
