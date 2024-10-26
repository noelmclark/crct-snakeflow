#### These rules use Plink2.0 & Plink1.9 on our filtered BCF ####

### PLINK based rules ###
## these rules take our filtered BCF and run PLINK2 on it to generate a PCA
# --allow-extra-chromosomes lets the bed file include non-human chroms (PLINK is defaulted to humans)
# I removed the Y chrom, mitogenome, and unassembled scaffolds using bcftools, but you could use the --not-chr flag in PLINK
# the --pheno option splits our bcf into one per population that includes all indivdiuals with the population id from the popfile
# this is so --freq can generate population specific allele frequencies vs ones based on the all sample bcf
rule calc_allele_freq:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf.csi",
        popfile="config/plink-popfile.tsv",
    output:
        afreq="results/plink/allele-freq/aut-snps-{maf}",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/allele-freq/aut-snps-{maf}.log",
    benchmark:
        "results/benchmarks/plink/allele-freq/aut-snps-{maf}.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --pheno {input.popfile} "
        " --allow-extra-chr "
        " --freq --loop-cats population "
        " --out {output.afreq} 2> {log} "

## This rules generates a PCA using Plink2.0 from our filtered BCF
# the --geno 0.01 applies a 10% missingness threshold filter
rule make_plink_pca:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf.csi",
        afreq="results/plink/allele-freq/aut-snps-{maf}.afreq"
    output:
        pca="results/plink/pca/aut-snps-{maf}-pca",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/pca/aut-snps-{maf}-pca.log",
    benchmark:
        "results/benchmarks/plink/pca/aut-snps-{maf}-pca.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --read-freq {input.afreq} --set-missing-var-ids @:#[b37]\$r,\$a "
        " --pca --allow-extra-chr --geno 0.1 "
        " --out {output.pca} 2> {log} "


## This rule calculate pairwise Fst values using the Weir & Cockerham (1984) method 
# on our hard filtered BCF file with only biallelic snps that pass a MAF cutoff of 0.05
rule make_pw_fst_snp:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf",
        tbi="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf.csi",
        popfile="config/plink-popfile.tsv",
    output:
        fst="results/plink/pw-fst/aut-snps-{maf}-fst",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/pw-fst/aut-snps-{maf}-fst.log",
    benchmark:
        "results/benchmarks/plink/pw-fst/aut-snps-{maf}-fst.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --const-fid 0 --within {input.popfile} population"
        " --out {output.fst} --fst population method=wc 2> {log} "

## This rule generates a Phylip file from the filtered BCF which is used as input for the IQ Tree and splits-tree programs
# need to update this
rule make_phylip:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf",
        csi="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf.csi",
        popfile="config/plink-popfile.tsv",
    output:
        phylip="results/plink/phylip/all-samp-no-y",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/phylip/all-samp-no-y.log",
    benchmark:
        "results/benchmarks/plink/phylip/all-samp-no-y.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr --not-chr NC_048593.1 "
        " --geno 0.25 --snps-only "
        " --export phylip used-sites --out {output.phylip} 2> {log} "


### BONEYARD (i.e. rules not used currently)

## This rule identifies runs of homozygosity (ROH) using the PLINK1.9 method
# https://www.cog-genomics.org/plink/1.9/ibd#homozyg

## This rule calculate pairwise Fst values using the Weir & Cockerham (1984) method 
# on our hard filtered BCF file with all variants
# populations are loaded as categorical phenotypes first from a popfile
# if we want to do jacknife blocksize as well: https://groups.google.com/g/plink2-users/c/vr8rHzYVhZo/m/syHgeWLYAAAJ
# I don't call for this currently 
rule make_pw_fst_all:
    input:
        bcf="results/bcf/all.bcf",
        tbi="results/bcf/all.bcf.csi",
        popfile="config/plink-popfile.tsv",
    output:
        fst="results/plink/pw-fst/fst-all",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=18800,
        cpus=2,
    threads: 2
    log:
        "results/logs/plink/pw-fst/fst-all.log",
    benchmark:
        "results/benchmarks/plink/pw-fst/fst-all.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --new-id-max-allele-len 387 "
        " --allow-extra-chr "
        " --const-fid 0 --within {input.popfile} population"
        " --out {output.fst} --fst population method=wc 2> {log} "