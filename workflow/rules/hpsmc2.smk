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
        " bcftools mpileup --full-BAQ -Ou -f {input.ref} -q30 -Q60 -r {wildcards.hpsmcchroms} {input.bam} | "
        " {input.dir}/Chrom-Compare/pu2fa -c {wildcards.hpsmcchroms} -C 50 -b -s 1 -e {params.end} > {output} 2> {log} "


rule concat_haploidized_bam:
    input:
        expand("results/hpsmc-test/haploidize_bam_sect/{{hpsmcpops}}/{c}_haploidized.fa", c=hpsmcchroms),
    output:
        "results/hpsmc-test/haploidized_bam/{hpsmcpops}_haploidized.fa",
    log:
        "results/logs/hpsmc-test/concat_haploidized_bam/{hpsmcpops}.log",
    benchmark:
        "results/benchmarks/hpsmc-test/concat_haploidized_bam/{hpsmcpops}.log",
    shell:
        " cat {input} > {output} 2> {log} "


rule psmcfa_from_2_fastas:
    input:
        pop1="results/hpsmc-test/haploidized_bam/{pop1}_haploidized.fa",
        pop2="results/hpsmc-test/haploidized_bam/{pop2}_haploidized.fa"
    output:
        "results/hpsmc-test/psmcfa-from-2-fastas/{pop1}---x---{pop2}.psmcfa"
    conda:
        "../envs/hpsmc.yaml"
    resources:
        time="23:59:59",
        mem_mb=9400,
        cpus=2,
    log:
        "results/logs/hpsmc-test/psmcfa-from-2-fastas/{pop1}---x---{pop2}.log"
    benchmark:
        "results/benchmarks/hpsmc-test/psmcfa-from-2-fastas/{pop1}---x---{pop2}.bmk"
    shell:
        "python workflow/scripts/hPSMC/JO_psmcfa_from_2_fastas.py -b10 -m5 {input.pop1} {input.pop2} > {output} 2> {log}"