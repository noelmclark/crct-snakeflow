#### This rule file is used to subset our original BCF (aut-bisnp-no5indel.bcf) to just
### samples from the SRM and their associated variants - identified as necessary following defense

# this rule generates a global allele freq file
# that is necessary for running PCA with all samples 
rule make_srm_pgen:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi",
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
        " --remove-if population == coastal "
        " --remove-if population == westslope "
        " --remove-if population == lahontan "
        " --remove-if population == lahontan "
        " --make-pgen "
        " --out {output.pgen} 2> {log} "
