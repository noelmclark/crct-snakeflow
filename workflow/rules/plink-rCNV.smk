#### These rules use Plink2.0 & Plink1.9 on our BCFs from the rCNV method ####
## --allow-extra-chromosomes lets the bed file include non-human chroms (PLINK is defaulted to humans)
## with aut-bcf files, I removed the Y chrom, mitogenome, and unassembled scaffolds using bcftools,
## but you could use the --not-chr flag in PLINK##


#########################################################
## this chunk uses aut-bisnps-no5indel-rcnv-by-dvs-2.0.bcf ##
#########################################################

# this rule generates a global allele freq file that is necessary for running PCA with all samples 
rule calc_allele_counts_rcnv_dvs:
    input:
        bcf="results/bcf/aut-bisnps-no5indel-rcnv-by-dvs-2.0.bcf",
        tbi="results/bcf/aut-bisnps-no5indel-rcnv-by-dvs-2.0.bcf.csi",
    output:
        acount="results/plink/allele-count/aut-bisnps-no5indel-rcnv-by-dvs-2.0",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=7480    
    log:
        "results/logs/plink/allele-count/aut-bisnps-no5indel-rcnv-by-dvs-2.0.log",
    benchmark:
        "results/benchmarks/plink/allele-count/aut-bisnps-no5indel-rcnv-by-dvs-2.0.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --freq counts "
        " --out {output.acount} 2> {log} "



## this rule removes variants that don't pass a 10% missingness filter
# and creates a new PLINK2 binary fileset (pgen) that can be used for downstream rules
rule filter_missingness_rcnv_dvs:
    input:
        bcf="results/bcf/aut-bisnps-no5indel-rcnv-by-dvs-2.0.bcf",
        tbi="results/bcf/aut-bisnps-no5indel-rcnv-by-dvs-2.0.bcf.csi",
    output:
        "results/plink/missingness/aut-bisnps-no5indel-rcnv-by-dvs-2.0",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=7480 
    log:
        "results/logs/plink/missingness/aut-bisnps-no5indel-rcnv-by-dvs-2.0.log",
    benchmark:
        "results/benchmarks/plink/missingness/aut-bisnps-no5indel-rcnv-by-dvs-2.0.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --make-pgen "
        " --out {output} 2> {log} "


## This rules generates a .bed, .bim, and .fam file (input needed for running ADMIXTURE) 
# that includes only sites that pass our 10% missingness filter 
rule make_plink_bed_rcnv_dvs:
    input:
        acount="results/plink/allele-count/aut-bisnps-no5indel-rcnv-by-dvs-2.0.acount",
        popfile="config/plink-popfile.tsv",
        outlier="config/plink-outlier.tsv",
    params:
        pfile="results/plink/missingness/aut-bisnps-no5indel-rcnv-by-dvs-2.0",
    output:
        bed="results/plink/bed/MAC5/aut-bisnps-no5indel-rcnv-by-dvs-2.0-nooutlier-MAC5",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=7480
    log:
        "results/logs/plink/bed/MAC5/aut-bisnps-no5indel-rcnv-by-dvs-2.0-nooutlier-MAC5.log",
    benchmark:
        "results/benchmarks/plink/bed/MAC5/aut-bisnps-no5indel-rcnv-by-dvs-2.0-nooutlier-MAC5.bmk",
    shell:
        " plink2 --pfile {params.pfile} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --remove {input.outlier} "
        " --read-freq {input.acount} "
        " --mac 6 "
        " --pheno {input.popfile} "
        " --make-bed "
        " --out {output.bed} 2> {log} "



## This rules generates a PCA from our filtered BCF
# the --geno 0.01 applies a 10% missingness threshold filter  
rule make_plink_pca_rcnv_dvs:
    input:
        acount="results/plink/allele-count/aut-bisnps-no5indel-rcnv-by-dvs-2.0.acount",
        outlier="config/plink-outlier.tsv",
    params:
        pfile="results/plink/missingness/aut-bisnps-no5indel-rcnv-by-dvs-2.0",
    output:
        pca="results/plink/pca/MAC5/aut-bisnps-no5indel-rcnv-by-dvs-2.0-nooutlier-MAC5-pca",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=11220
    log:
        "results/logs/plink/pca/MAC5/aut-bisnps-no5indel-rcnv-by-dvs-2.0-nooutlier-MAC5-pca.log",
    benchmark:
        "results/benchmarks/plink/pca/MAC5/aut-bisnps-no5indel-rcnv-by-dvs-2.0-nooutlier-MAC5-pca.bmk",
    shell:
        " plink2 --pfile {params.pfile} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --remove {input.outlier} "
        " --read-freq {input.acount} "
        " --mac 6 "
        " --pca "
        " --out {output.pca} 2> {log} "



## This rule generates a Phylip file from the filtered BCF which is used as input for the IQ Tree and splits-tree programs
# the --geno 0.01 applies a 10% missingness threshold filter that should be redundant when using the new purned sites 
# need to double check file options and recommended filters in PLINK and IQTree
rule make_phylip_rcnv_dvs:
    input:
        acount="results/plink/allele-count/aut-bisnps-no5indel-rcnv-by-dvs-2.0.acount",
        outlier="config/plink-outlier.tsv",
    params:
        pfile="results/plink/missingness/aut-bisnps-no5indel-rcnv-by-dvs-2.0",
    output:
        phylip="results/plink/phylip/MAC5/aut-bisnps-no5indel-rcnv-by-dvs-2.0-nooutlier-MAC5",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=11220
    log:
        "results/logs/plink/phylip/MAC5/aut-bisnps-no5indel-rcnv-by-dvs-2.0-nooutlier-MAC5.log",
    benchmark:
        "results/benchmarks/plink/phylip/MAC5/aut-bisnps-no5indel-rcnv-by-dvs-2.0-nooutlier-MAC5.bmk",
    shell:
        " plink2 --pfile {params.pfile} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --remove {input.outlier} "
        " --read-freq {input.acount} "
        " --mac 6 "
        " --snps-only "
        " --export phylip used-sites "
        " --out {output.phylip} 2> {log} "


## this rule breaks the original bcf input file into individuals using --indv
# and produces a genotype count report that I can extract #het/#calls from
# we don't filter on missingess here because we're doing it by individual 
# eventually we will filter on individual level MACs
rule make_gt_count_rcnv_dvs:
    input:
        bcf="results/bcf/aut-bisnps-no5indel-rcnv-by-dvs-2.0.bcf",
        tbi="results/bcf/aut-bisnps-no5indel-rcnv-by-dvs-2.0.bcf.csi",
        acount="results/plink/allele-count/aut-bisnps-no5indel-rcnv-by-dvs-2.0.acount",
    output:
        gcount="results/plink/gt-count/MAC5/{sample}-aut-bisnps-no5indel-rcnv-by-dvs-2.0-MAC5",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=11220
    log:
        "results/logs/plink/gt-count/MAC5/{sample}-aut-bisnps-no5indel-rcnv-by-dvs-2.0-MAC5.log",
    benchmark:
        "results/benchmarks/plink/gt-count/MAC5/{sample}-aut-bisnps-no5indel-rcnv-by-dvs-2.0-MAC5.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --indv {wildcards.sample} "
        " --read-freq {input.acount} "
        " --mac 6 "
        " --geno-counts "
        " --out {output.gcount} 2> {log} "



###############################################################################



#########################################################
## this chunk uses aut-bisnps-no5indel-rcnv-by-cnv-2.0.bcf ##
#########################################################

# this rule generates a global allele freq file that is necessary for running PCA with all samples 
rule calc_allele_counts_rcnv_cnv:
    input:
        bcf="results/bcf/aut-bisnps-no5indel-rcnv-by-cnv-2.0.bcf",
        tbi="results/bcf/aut-bisnps-no5indel-rcnv-by-cnv-2.0.bcf.csi",
    output:
        acount="results/plink/allele-count/aut-bisnps-no5indel-rcnv-by-cnv-2.0",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=7480    
    log:
        "results/logs/plink/allele-count/aut-bisnps-no5indel-rcnv-by-cnv-2.0.log",
    benchmark:
        "results/benchmarks/plink/allele-count/aut-bisnps-no5indel-rcnv-by-cnv-2.0.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --freq counts "
        " --out {output.acount} 2> {log} "



## this rule removes variants that don't pass a 10% missingness filter
# and creates a new PLINK2 binary fileset (pgen) that can be used for downstream rules
rule filter_missingness_rcnv_cnv:
    input:
        bcf="results/bcf/aut-bisnps-no5indel-rcnv-by-cnv-2.0.bcf",
        tbi="results/bcf/aut-bisnps-no5indel-rcnv-by-cnv-2.0.bcf.csi",
    output:
        "results/plink/missingness/aut-bisnps-no5indel-rcnv-by-cnv-2.0",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=7480 
    log:
        "results/logs/plink/missingness/aut-bisnps-no5indel-rcnv-by-cnv-2.0.log",
    benchmark:
        "results/benchmarks/plink/missingness/aut-bisnps-no5indel-rcnv-by-cnv-2.0.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --make-pgen "
        " --out {output} 2> {log} "


## This rules generates a .bed, .bim, and .fam file (input needed for running ADMIXTURE) 
# that includes only sites that pass our 10% missingness filter 
rule make_plink_bed_rcnv_cnv:
    input:
        acount="results/plink/allele-count/aut-bisnps-no5indel-rcnv-by-cnv-2.0.acount",
        popfile="config/plink-popfile.tsv",
        outlier="config/plink-outlier.tsv",
    params:
        pfile="results/plink/missingness/aut-bisnps-no5indel-rcnv-by-cnv-2.0",
    output:
        bed="results/plink/bed/MAC5/aut-bisnps-no5indel-rcnv-by-cnv-2.0-nooutlier-MAC5",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=7480
    log:
        "results/logs/plink/bed/MAC5/aut-bisnps-no5indel-rcnv-by-cnv-2.0-nooutlier-MAC5.log",
    benchmark:
        "results/benchmarks/plink/bed/MAC5/aut-bisnps-no5indel-rcnv-by-cnv-2.0-nooutlier-MAC5.bmk",
    shell:
        " plink2 --pfile {params.pfile} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --remove {input.outlier} "
        " --read-freq {input.acount} "
        " --mac 6 "
        " --pheno {input.popfile} "
        " --make-bed "
        " --out {output.bed} 2> {log} "



## This rules generates a PCA from our filtered BCF
# the --geno 0.01 applies a 10% missingness threshold filter  
rule make_plink_pca_rcnv_cnv:
    input:
        acount="results/plink/allele-count/aut-bisnps-no5indel-rcnv-by-cnv-2.0.acount",
        outlier="config/plink-outlier.tsv",
    params:
        pfile="results/plink/missingness/aut-bisnps-no5indel-rcnv-by-cnv-2.0",
    output:
        pca="results/plink/pca/MAC5/aut-bisnps-no5indel-rcnv-by-cnv-2.0-nooutlier-MAC5-pca",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=11220
    log:
        "results/logs/plink/pca/MAC5/aut-bisnps-no5indel-rcnv-by-cnv-2.0-nooutlier-MAC5-pca.log",
    benchmark:
        "results/benchmarks/plink/pca/MAC5/aut-bisnps-no5indel-rcnv-by-cnv-2.0-nooutlier-MAC5-pca.bmk",
    shell:
        " plink2 --pfile {params.pfile} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --remove {input.outlier} "
        " --read-freq {input.acount} "
        " --mac 6 "
        " --pca "
        " --out {output.pca} 2> {log} "



## This rule generates a Phylip file from the filtered BCF which is used as input for the IQ Tree and splits-tree programs
# the --geno 0.01 applies a 10% missingness threshold filter that should be redundant when using the new purned sites 
# need to double check file options and recommended filters in PLINK and IQTree
rule make_phylip_rcnv_cnv:
    input:
        acount="results/plink/allele-count/aut-bisnps-no5indel-rcnv-by-cnv-2.0.acount",
        outlier="config/plink-outlier.tsv",
    params:
        pfile="results/plink/missingness/aut-bisnps-no5indel-rcnv-by-cnv-2.0",
    output:
        phylip="results/plink/phylip/MAC5/aut-bisnps-no5indel-rcnv-by-cnv-2.0-nooutlier-MAC5",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=11220
    log:
        "results/logs/plink/phylip/MAC5/aut-bisnps-no5indel-rcnv-by-cnv-2.0-nooutlier-MAC5.log",
    benchmark:
        "results/benchmarks/plink/phylip/MAC5/aut-bisnps-no5indel-rcnv-by-cnv-2.0-nooutlier-MAC5.bmk",
    shell:
        " plink2 --pfile {params.pfile} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --remove {input.outlier} "
        " --read-freq {input.acount} "
        " --mac 6 "
        " --snps-only "
        " --export phylip used-sites "
        " --out {output.phylip} 2> {log} "


## this rule breaks the original bcf input file into individuals using --indv
# and produces a genotype count report that I can extract #het/#calls from
# we don't filter on missingess here because we're doing it by individual 
# eventually we will filter on individual level MACs
rule make_gt_count_rcnv_cnv:
    input:
        bcf="results/bcf/aut-bisnps-no5indel-rcnv-by-cnv-2.0.bcf",
        tbi="results/bcf/aut-bisnps-no5indel-rcnv-by-cnv-2.0.bcf.csi",
        acount="results/plink/allele-count/aut-bisnps-no5indel-rcnv-by-cnv-2.0.acount",
    output:
        gcount="results/plink/gt-count/MAC5/{sample}-aut-bisnps-no5indel-rcnv-by-cnv-2.0-MAC5",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=11220
    log:
        "results/logs/plink/gt-count/MAC5/{sample}-aut-bisnps-no5indel-rcnv-by-cnv-2.0-MAC5.log",
    benchmark:
        "results/benchmarks/plink/gt-count/MAC5/{sample}-aut-bisnps-no5indel-rcnv-by-cnv-2.0-MAC5.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --indv {wildcards.sample} "
        " --read-freq {input.acount} "
        " --mac 6 "
        " --geno-counts "
        " --out {output.gcount} 2> {log} "

###############################################################################