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

rule srm_allele_counts:
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

## This rules generates a .bed, .bim, and .fam file (input needed for running ADMIXTURE) 
# that includes only sites that pass our 10% missingness filter 
rule srm_plink_bed:
    input:
        acount="results/plink/srm-subset/allele-count/srm-aut-bisnps-no5indel.acount",
        popfile="config/plink-popfile.tsv",
        srmfile="config/plink-srm-samples.tsv"
    params:
        pfile="results/plink/srm-subset/make-pgen/srm-aut-bisnps-no5indel",
    output:
        bed="results/plink/srm-subset/bed/MAC1/srm-aut-bisnps-no5indel-MAC1",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=7480
    log:
        "results/logs/plink/srm-subset/bed/MAC1/srm-aut-bisnps-no5indel-MAC1.log",
    benchmark:
        "results/benchmarks/plink/srm-subset/bed/MAC1/srm-aut-bisnps-no5indel-MAC1.bmk",
    shell:
        " plink2 --pfile {params.pfile} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --read-freq {input.acount} "
        " --mac 2 "
        " --pheno {input.popfile} "
        " --make-bed "
        " --out {output.bed} 2> {log} "



## This rules generates a PCA from our filtered BCF
# the --geno 0.01 applies a 10% missingness threshold filter  
rule srm_plink_pca:
    input:
        acount="results/plink/srm-subset/allele-count/srm-aut-bisnps-no5indel.acount",
    params:
        pfile="results/plink/srm-subset/make-pgen/srm-aut-bisnps-no5indel",
    output:
        pca="results/plink/srm-subset/pca/MAC1/srm-aut-bisnps-no5indel-MAC1-pca",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=11220
    log:
        "results/logs/plink/srm-subset/pca/MAC1/srm-aut-bisnps-no5indel-MAC1-pca.log",
    benchmark:
        "results/benchmarks/plink/srm-subset/pca/MAC1/srm-aut-bisnps-no5indel-MAC1-pca.bmk",
    shell:
        " plink2 --pfile {params.pfile} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --read-freq {input.acount} "
        " --mac 2 "
        " --pca "
        " --out {output.pca} 2> {log} "