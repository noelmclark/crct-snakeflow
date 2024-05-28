#rule to calculate avg depth of coverage of bam file using samtools
rule get_coverage_depth:
    input:
        bam="results/mkdup/{sample}.bam"
    output:
        "results/qc/coverage/{sample}.txt"
    conda:
        "../envs/sambcftools.yaml"
    log:
        "results/get-coverage-depth/{sample}.log"
    benchmark:
        "results/benchmarks/get-coverage-depth/{sample}.bmk"
    shell:
        "samtools depth -a -H {input.bam} -o {output} 2> {log}"
