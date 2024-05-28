rule vcf_stats:
    input:
        "results/vcf/all.vcf.gz"
    output:
        "results/vcf-stats/all.vcf.stats"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/vcf-stats/all.log"
    benchmark:
        "results/benchmarks/all-vcf-stats.bmk"
    shell:
        "bcftools stats {input} > {output} 2> {log}"