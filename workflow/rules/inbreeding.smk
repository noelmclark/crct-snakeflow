## this rule counts the total variable sites needed before counting hets in next rule
rule bcftools_stats:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf",
        csi="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf.csi",
    output:
        stats="results/inbreeding/het/{sample}-stats.txt",
    conda:
        "./envs/bcftools.yaml",
    log:
        "results/inbreeding/het/{sample}-stats.log",
    benchmark:
        "results/inbreeding/het/{sample}-stats.bmk"
    shell:
        " bcftools stats -s {wildcards.sample} {input} > {output.stat} 2> {log} "

## This rule counts the heterozygous sites (0/1) in an individual vcf
rule count_hets:
    input:
        "results/inbreeding/het/{sample}-stats.txt",
    output:
        "results/inbreeding/het/{sample}-het.txt",
    log:
        "results/inbreeding/het/{sample}-het.log",
    benchmark:
        "results/inbreeding/het/{sample}-het.bmk"
    shell:
        " grep \"0/1:\" {input} | wc -l > {output} 2> {log} "




## Not currently used

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
        tbi="results/roh/allele-freq/snps-maf-0.05/{population}-freqs.tabs.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    log:
        freq="results/logs/roh/allele-freq/snps-maf-0.05/{population}-freqs.log",
        index="results/logs/roh/allele-freq/snps-maf-0.05/{population}-freqs.tbi.log"
    benchmark:
        "results/benchmarks/roh/allele-freq/snps-maf-0.05/{population}-freqs.bmk",
    shell:
        " bcftools view -s {params.pops} {input.bcf} |"
        " bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' | "
        " bgzip -c > {output.afreq} 2> {log.freq} && "
        " tabix -f -s1 -b2 -e2 {output.afreq} > {output.tbi} 2> {log.index} "


rule run_bcftools_roh:
    input:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
        afreq="results/roh/allele-freq/snps-maf-0.05/{population}-freqs.tabs.gz",
        tbi="results/roh/allele-freq/snps-maf-0.05/{population}-freqs.tabs.gz.tbi"
    params:
        pops=get_comma_sep_pop_names
    output:
        "results/roh/snps-maf-0.05/{population}-roh.txt",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/roh/snps-maf-0.05/{population}-roh.log",
    benchmark:
        "results/benchmarks/roh/snps-maf-0.05/{population}-roh.bmk",
    shell:
        " bcftools view -s {params.pops} {input.bcf} | "
        " bcftools roh --AF-file {input.afreq} --GTs-only 30 "
        " -o {output} 2> {log} "