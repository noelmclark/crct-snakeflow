rule test_haploidize_bam_sect:
    input:
        bam=get_hpsmc_bams_in_pop,
        ref="resources/genome/OmykA.fasta",
        dir="results/chromcompare"
    output:
        "results/hpsmc-test/haploidize_bam_sect/{hpsmcpops}/{hpsmcchroms}_haploidized.fa",
    params:
        end=get_hpsmc_chrom_end
    conda:
        "../envs/bcftools-chromcompare.yaml"
    log:
        "results/logs/hpsmc-test/haploidize-bam-sect/{hpsmcpops}/{hpsmcchroms}.log",
    benchmark:
        "results/benchmarks/hpsmc-test/haploidize-bam-sect/{hpsmcpops}/{hpsmcchroms}.bmk",
    shell:
        " bcftools mpileup --full-BAQ -Ou -f {input.ref} -r {wildcards.hpsmcchroms} {input.bam} | "
        " {input.dir}/Chrom-Compare/pu2fa -c {wildcards.hpsmcchroms} -C 50 -s 1 -e {params.end} > {output} 2> {log} "