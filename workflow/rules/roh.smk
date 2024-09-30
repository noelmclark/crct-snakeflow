# this filters the bcf file by sample names listed in the given population wildcard 
# and generates an allele frequency file for each population
rule make_bcftools_pop_afreq:
    input:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
    params:
        pops=get_comma_sep_pop_names
    output:
        afreq="results/roh/allele-freq/snps-maf-0.05/{population}-freqs.tabs.gz",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/roh/allele-freq/snps-maf-0.05/{population}-freqs.log",
    benchmark:
        "results/benchmarks/roh/allele-freq/snps-maf-0.05/{population}-freqs.bmk",
    shell:
        " bcftools view -s {params.pops} {input.bcf} | "
        " bctools query -f'%CHROM\t%POS\t%REF,%ALT\t[%AF]\n' | "
        " bgzip -c > {output.afreq} 2> {log} "


rule run_bcftools_roh:
    input:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
        afreq="results/roh/allele-freq/snps-maf-0.05/{population}-freqs.tabs.gz",
    output:
        "results/roh/snps-maf-0.05-roh.?",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/roh/snps-maf-0.05-roh.log",
    benchmark:
        "results/benchmarks/roh/snps-maf-0.05-roh.bmk",
    shell:
        " bcftools view -s {params.pops} {input.bcf} | "
        " bcftools roh --AF-file {input.afreq} --GTs-only 30 "
        " -o {output} 2> {log} "