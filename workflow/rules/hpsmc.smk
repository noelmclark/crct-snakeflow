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

rule test_haploidize_bam:
    input:
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",
        ref="resources/genome/OmykA.fasta",
        sgc={sg_or_chrom}
    output:
        "results/hpsmc/test_haploidize_bam/{sample}_test_haploid.txt"
    log:
        "results/logs/hpsmc/test_haploidize_bam/{sample}.log",
    shell:
        " for i in {input.sgc}; do "
        "  echo 'bcftools mpileup --full-BAQ -s -Ou -f {input.ref} -q30 -Q60 -r $i {input.bam} | "
        "  pu2fa -c $i -C 50' ; done > {output} 2> {log} "
