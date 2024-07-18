### These rules generate a PCA using Plink2.0 from our filtered BCF

## The following 2 rules are copied from Eric's post-bcf workflow (bcftools_filter.smk) with edits
# https://github.com/eriqande/mega-post-bcf-exploratory-snakeflows

rule bcf_filt_scatter:
    input:
        bcf="results/bcf/all.bcf",
        csi="results/bcf/all.bcf.csi",
        sfile="config/sample_subsets/all-fish.txt"
        scat_path="results/scatter_config/scatters_1200000.tsv"
    params:
        bcftools_opts="-v snps -m 2 -M 2 -i 'FILTER=\"PASS\" && MAF >= 0.05'",
        scat="{scatter}"
    output:
        scat_regions="results/pca/scat_regions/{scatter}.scat_regions.tsv",
        bcf="results/bcf/filt_biallelic_maf_0.05/sections/{scatter}.bcf",
        stats="results/bcf/filt_biallelic_maf_0.05/sections/{scatter}.bcf_stats.txt",
        pos="results/bcf/filt_biallelic_maf_0.05/sections/{scatter}.positions.tsv.gz",
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
        "( bcftools concat --naive {input.bcfs} > {output.bcf} && "
        "  bcftools index {output.bcf} && "
        "  cat {input.poses} > {output.pos} && "
        "  plot-vcfstats -m {input.statses} > {output.stats} "
        ") 2> {log} "