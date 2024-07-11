### These rules implement the hybrid PSMC method explained in Cahill et al. 2016 (https://doi.org/10.1098/rstb.2015.0138)
## and documented in https://github.com/jacahill/hPSMC
## 1. Create an hPSMC.psmcfa input file from two samples
## 2. run PSMC using the hPSMC.psmcfa
## 3. visualize hPSMC.psmc and estimate pre-divergence pop. size
## 4. Run simulations of divergence without post-divergence migration to compare to the hPSMC plot
## 5. Plot simulations with the original data to show divergence between samples
## Important Warning: By default ms outputs 5 decimal places for mutation locations which is enough to 100,000 bins. 
# Recent versions of ms include the -p flag which allows you to set the number of decimal places to report. 
# I recommend using -p8 in most cases.

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


rule test_haploidize_bam:
    input:
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",
        bai="results/angsd_bams/overlap_clipped/{sample}.bai",
        ref="resources/genome/OmykA.fasta",
        idx="resources/genome/OmykA.dict",
        fai="resources/genome/OmykA.fasta.fai",
        sgc={sg_or_chrom}
    output:
        "results/hpsmc/test_haploidize_bam/{sample}_test_haploid.txt",
    log:
        "results/logs/hpsmc/test_haploidize_bam/{sample}/{sg_or_chrom}.log",
    shell:
        " for i in {input.sgc}; do "
        "  echo 'bcftools mpileup --full-BAQ -s -Ou -f {input.ref} -q30 -Q60 -r $i {input.bam} | "
        "  pu2fa -c $i -C 50' ; done > {output} 2> {log} "

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