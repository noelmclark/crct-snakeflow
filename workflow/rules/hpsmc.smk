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

## rule to install chromcompare into the bcftools-chromcompare.yaml
rule install_chromcompare:
    params:
        url=config["params"]["chromcompare"]["url"]
    output:  
        flagfile=touch("results/flags/chromcompare_installed")
    conda:
        "../envs/bcftools-chromcompare.yaml"
    log:
        "results/logs/install_chromcompare/log.txt"
    shell:
        "(TMP=$(mktemp -d) && cd $TMP && "
        " git clone {params.url} && "
        " cd Chrom-Compare && "
        " make ) > {log} 2>&1  "

rule haploidize_bam_sections:
    input:
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",
        ref="resources/genome/OmykA.fasta",
    output:
        temp("results/hpsmc/haploidize_bam_sect/{sample}/{chromo}_haploidized.fa"),
    conda:
        "../envs/bcftools-chromcompare.yaml"
    resources:
        time="23:59:59",
    log:
        "results/logs/hpsmc/haploidize-bam-sect/{sample}/{chromo}.log",
    benchmark:
        "results/benchmarks/hpsmc/haploidize-bam-sect/{sample}/{chromo}.bmk",
    shell:
        " echo 'bcftools mpileup --full-BAQ -s -Ou -f {input.ref} -q30 -Q60 -r {wildcards.chromo} {input.bam} | "
        " pu2fa -c {wildcards.chromo} -C 50 > {output}' "

rule concat_haploidized_bam:
    input:
        expand("results/hpsmc/haploidize_bam_sect/{{sample}}/{c}_haploidized.fa", c=unique_chromosomes),
    output:
        "results/hpsmc/haploidized_bam/{sample}_haploidized.fa",
    log:
        "results/logs/hpsmc/concat_haploidized_bam/{sample}.log",
    benchmark:
        "results/benchmarks/hpsmc/concat_haploidized_bam/{sample}.log",
    shell:
        " cat {input} > {output} 2> {log} "