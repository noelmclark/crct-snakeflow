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


#rule to run stamtools stats on each of the final mkdup bams
rule samtools_stats:
    input:
        "results/mkdup/{sample}.bam",
    output:
        "results/qc/samtools_stats/{sample}.txt",
    log:
        "results/logs/samtools_stats/{sample}.log",
    benchmark:
        "results/benchmarks/samtools_stats/{sample}.bmk",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input} > {output} 2> {log} "