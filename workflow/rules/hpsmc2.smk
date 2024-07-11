rule haploidize_bam_sections:
    input:
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",
        bai="results/angsd_bams/overlap_clipped/{sample}.bai",
        ref="resources/genome/OmykA.fasta",
        idx="resources/genome/OmykA.dict",
        fai="resources/genome/OmykA.fasta.fai",
        sgc={sg_or_chrom}
    output:
        "results/hpsmc/haploidize_bam_sect/{sample}/{sg_or_chrom}_haploidized.fa",
    conda:
        "../envs/bcftools-pu2fa.yaml"
    #resources:
    log:
        "results/logs/hpsmc/haploidize-bam-sect/{sample}/{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hpsmc/haploidize-bam-sect/{sample}/{sg_or_chrom}.bmk",
    shell:
        " bcftools mpileup --full-BAQ -s -Ou -f {input.ref} -q30 -Q60 -r {input.sgc} {input.bam} | "
        " pu2fa -c {input.sgc} -C 50 > {output} "

rule concat_haploidized_bam:
    input:
        expand("results/hpsmc/haploidize_bam_sect/{{sample}}/{sgc}_haploidized.fa", sgc=sg_or_chrom),
    output:
        "results/hpsmc/haploidized_bam/{sample}_haploidized.fa",
    log:
        "results/logs/hpsmc/concat_haploidized_bam/{sample}.log",
    benchmark:
        "results/benchmarks/hpsmc/concat_haploidized_bam/{sample}.log",
    shell:
        " cat {input} > {output} 2> {log} "


#is there a way to do both rules at once?
#        " for i in {input.sgc}; do "
#        "  bcftools mpileup --full-BAQ -s -Ou -f {input.ref} -q30 -Q60 -r $i {input.bam} | "
#        "  pu2fa -c $i -C 50 >> {output} 2>> {log}; "
#        " done "
#
#        or
#
#
#       " for i in {input.sgc}; do "
#        "  bcftools mpileup --full-BAQ -s -Ou -f {input.ref} -q30 -Q60 -r $i {input.bam} | "
#        "  pu2fa -c $i -C 50 ; done > {output} 2> {log} "