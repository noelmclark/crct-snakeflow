#### These rules use Plink2.0 & Plink1.9 on our filtered BCF ####

### The following 2 rules are copied from Eric's post-bcf workflow (bcftools_filter.smk) with edits
## https://github.com/eriqande/mega-post-bcf-exploratory-snakeflows
## they generate a BCF file with only biallelic snps that pass a MAF cuttoff of 0.05
rule bcf_filt_scatter:
    input:
        bcf="results/bcf/all.bcf",
        csi="results/bcf/all.bcf.csi",
        sfile="config/sample_subsets/all-fish.txt",
        scat_path="results/scatter_config/scatters_1200000.tsv"
    params:
        bcftools_opts="-v snps -m 2 -M 2 -i 'FILTER=\"PASS\" && MAF >= 0.05'",
        scat="{scatter}"
    output:
        scat_regions=temp("results/pca/scat_regions/{scatter}.scat_regions.tsv"),
        bcf=temp("results/bcf/filt_biallelic_maf_0.05/sections/{scatter}.bcf"),
        stats=temp("results/bcf/filt_biallelic_maf_0.05/sections/{scatter}.bcf_stats.txt"),
        pos=temp("results/bcf/filt_biallelic_maf_0.05/sections/{scatter}.positions.tsv.gz"),
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/bcf_filt_scatter/bcf/filt_biallelic_maf_0.05/sections/{scatter}.log"
    benchmark:
        "results/benchmarks/bcf_filt_scatter/filt_biallelic_maf_0.05/sections/{scatter}.bmk"
    shell:
        " (awk -v scat='{params.scat}' -f workflow/scripts/pca/get_scat_regions.awk {input.scat_path} > {output.scat_regions} && "
        "    bcftools view -Ou -R {output.scat_regions} {input.bcf} | "
        "    bcftools +fill-tags -Ou -- -t all | "
        "    bcftools view -Ob {params.bcftools_opts} > {output.bcf}  && "
        "    bcftools stats --af-tag MAF {output.bcf} > {output.stats} && "
        "    bcftools query -f '%CHROM\\t%POS\\n' {output.bcf} | gzip -c > {output.pos} "
        " ) 2> {log} "

rule bcf_filt_gather:
    input:
        bcfs=expand("results/bcf/filt_biallelic_maf_0.05/sections/{scat}.bcf", scat=unique_scats),
        poses=expand("results/bcf/filt_biallelic_maf_0.05/sections/{scat}.positions.tsv.gz", scat=unique_scats),
        statses=expand("results/bcf/filt_biallelic_maf_0.05/sections/{scat}.bcf_stats.txt", scat=unique_scats),
    output:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
        pos="results/bcf/filt_biallelic_maf_0.05/info/positions.tsv.gz",
        stats="results/bcf/filt_biallelic_maf_0.05/info/bcf_stats.txt",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/bcf/filt_biallelic_maf_0.05/main.log"
    benchmark:
        "results/benchmarks/bcf/filt_biallelic_maf_0.05/main.bmk"
    shell:
        " ( bcftools concat --naive {input.bcfs} | "
        "  bcftools sort -Ob > {output.bcf} && "
        "  bcftools index {output.bcf} && "
        "  cat {input.poses} > {output.pos} && "
        "  plot-vcfstats -m {input.statses} > {output.stats} "
        " ) 2> {log} "


### PLINK based rules ###
## these rules take our filtered BCF and run PLINK2 on it to generate a PCA
# --allow-extra-chromosomes lets the bed file include non-human chroms (PLINK is defaulted to humans)
# the --not-chr options removes the Y chromosome from the bed file
rule calc_allele_freq:
    input:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
    output:
        afreq="results/plink/allele-freq/snps-no-y",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/allele-freq/snps-no-y.log",
    benchmark:
        "results/benchmarks/plink/allele-freq/snps-no-y.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --freq --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr --not-chr NC_048593.1 "
        " --out {output.afreq} 2> {log} "

## This rules generates a PCA using Plink2.0 from our filtered BCF
rule make_plink_pca:
    input:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
        afreq="results/plink/allele-freq/snps-no-y.afreq"
    output:
        pca="results/plink/pca/snps-no-y-pca",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/pca/snps-no-y-pca.log",
    benchmark:
        "results/benchmarks/plink/pca/snps-no-y-pca.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --read-freq {input.afreq} --set-missing-var-ids @:#[b37]\$r,\$a "
        " --pca --allow-extra-chr --not-chr NC_048593.1 "
        " --out {output.pca} 2> {log} "

## This rule calculate pairwise Fst values using the Weir & Cockerham (1984) method 
# on our hard filtered BCF file with all variants
# populations are loaded as categorical phenotypes first from a popfile
# if we want to do jacknife blocksize as well: https://groups.google.com/g/plink2-users/c/vr8rHzYVhZo/m/syHgeWLYAAAJ
rule make_pw_fst_all:
    input:
        bcf="results/bcf/all.bcf",
        tbi="results/bcf/all.bcf.csi",
        popfile="config/plink-popfile.tsv",
    output:
        fst="results/plink/pw-fst-all/fst-all",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=18800,
        cpus=2,
    threads: 2
    log:
        "results/logs/plink/pw-fst-all/fst-all.log",
    benchmark:
        "results/benchmarks/plink/pw-fst-all/fst-all.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --new-id-max-allele-len 387 "
        " --allow-extra-chr --not-chr NC_048593.1 "
        " --const-fid 0 --within {input.popfile} population"
        " --out {output.fst} --fst population method=wc 2> {log} "

## This rule calculate pairwise Fst values using the Weir & Cockerham (1984) method 
# on our hard filtered BCF file with only biallelic snps that pass a MAF cutoff of 0.05
rule make_pw_fst_snp:
    input:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
        popfile="config/plink-popfile.tsv",
    output:
        fst="results/plink/pw-fst-snps-0.05/fst-snps-0.05",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/pw-fst-snps-0.05/fst-snps-0.05.log",
    benchmark:
        "results/benchmarks/plink/pw-fst-snps-0.05/fst-snps-0.05.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr --not-chr NC_048593.1 "
        " --const-fid 0 --within {input.popfile} population"
        " --out {output.fst} --fst population method=wc 2> {log} "

## This rule calculates a sample based missing data report from a BCF file

## This rule identifies runs of homozygosity (ROH) using the PLINK1.9 method
# https://www.cog-genomics.org/plink/1.9/ibd#homozyg