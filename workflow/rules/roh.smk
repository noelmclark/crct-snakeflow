# this is wrong right now because it assumes all individuals in the BCF are from the same population to generate an allele frequency file
# need to split up the bcf by population first...
make_bcftools_afreq:
    input:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
    output:
        afreq="results/roh/allele-freq/snps-maf-0.05-freqs.tabs.gz",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/roh/allele-freq/snps-maf-0.05-freqs.log",
    benchmark:
        "results/benchmarks/roh/allele-freq/snps-maf-0.05-freqs.bmk",
    shell:
        " bctools query -f'%CHROM\t%POS\t%REF,%ALT\t[%AF]\n' {input.bcf} "
        " bgzip -c > {output.afreq} 2> {log} "


run_bcftools_roh:
    input:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
        afreq="results/roh/allele-freq/snps-maf-0.05-freqs.tabs.gz",
    output:
        "results/roh/snps-maf-0.05-roh.?",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/roh/snps-maf-0.05-roh.log",
    benchmark:
        "results/benchmarks/roh/snps-maf-0.05-roh.bmk",
    shell:
        " bcftools roh --AF-file {input.afreq} --GTs-only 30 "
        " {input.bcf} -o {output}"