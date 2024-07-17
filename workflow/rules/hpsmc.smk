### These rules implement the hybrid PSMC method explained in Cahill et al. 2016 (https://doi.org/10.1098/rstb.2015.0138)
## and documented in https://github.com/jacahill/hPSMC
## 1. Create an hPSMC.psmcfa input file from two samples
## 2. run PSMC using the hPSMC.psmcfa
## 3. visualize hPSMC.psmc and estimate pre-divergence Ne & upper and lower divergence time estimates
## 4. Run simulations of divergence without post-divergence migration to compare to the original hPSMC plot
## 5. Plot simulations with the original data to show divergence between samples
## Important Warning: By default ms outputs 5 decimal places for mutation locations which is enough to 100,000 bins. 
# Recent versions of ms include the -p flag which allows you to set the number of decimal places to report. 
# I recommend using -p8 in most cases.

## this is a rule that will install Chrom-Compare from https://github.com/Paleogenomics/Chrom-Compare 
# into the active conda env and make (compile) it
# the pu2fa function from Chrom-Compare is needed to haploidize bam sections for hPMSC
# so, whenever you need to use Chrom-Compare, you call the same
# conda environment, and have the flagfile as an input
# depenency to ensure that Chrom-Compare has already been successfully
# built into that conda env. Similar to process for pcangsd from Eric's bcf snakeflow
rule install_chromcompare:
    params:
        hash=config["chromcompare"]["version"],
        url=config["chromcompare"]["url"]
    output:  
        flagfile=touch("results/flags/chromcompare_installed")
    conda:
        "../envs/bcftools-chromcompare.yaml"
    log:
        "results/logs/install_chromcompare/log.txt"
    shell:
        "(TMP=$(mktemp -d) && cd $TMP && "
        " git clone {params.url} && "
        " cd Chrom-Compare  && "
        " git checkout {params.hash} && "
        " make ) > {log} 2>&1  "


## 1. Create an hPSMC.psmcfa file for each combination of 2 samples 
rule haploidize_bam_sections:
    input:
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",
        ref="resources/genome/OmykA.fasta",
        flagfile="results/flags/chromcompare_installed"
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


#rule psmcfa_from_2_fastas:
#    input:
#        pop1="results/hpsmc/haploidized_bam/{pop1}_haploidized.fa",
#        pop2="sample2_all.fa"
#    output:
#        "results/hpsmc/psmcfa-from-2-fastas/{pop1}---x---{pop2}.hPSMC.psmcfa"
#    conda:
#        "../envs/hpsmc.yaml"
#    log:
#        ""
#    benchmark:
#        ""
#    shell:
#        "python workflow/scripts/hPSMC/psmcfa_from_2_fastas.py -b10 -m5 {input.pop1} {input.pop2} > {output}"


## 2. run PSMC on each of the pop1---x---pop2-hPSMC.psmcfa files
rule run_hpsmc:
    input:
        "results/hpsmc/psmcfa-from-2-fastas/{pop1}---x---{pop2}.hPSMC.psmcfa"
    output:
        "results/hpsmc/run-hpsmc/{pop1}---x---{pop2}.hPSMC.psmc"
    conda:
        "../envs/hpsmc.yaml"
    log:
        "results/logs/hpsmc/run-hpsmc/{pop1}---x---{pop2}.log"
    benchmark:
        "results/benchmarks/hpsmc/run-hpsmc/{pop1}---x---{pop2}.bmk"
    shell:
        "psmc -N25 -t15 -r5 -p '4+25*2+4+6' -o {output} {input} 2> {log}"

## 3. visualize hPSMC plots (using PSMC) and esimate pre-divergence Ne & upper and lower divergence time by looking at plots
## rule to plot hpsmc to visualize result
# -u [per-generation mutation rate] from https://doi.org/10.1371/journal.pgen.1010918
# -g [generation time in years] 
rule hpsmc_plot:
    input:
        "results/hpsmc/run-hpsmc/{pop1}---x---{pop2}.hPSMC.psmc"
    output:
        eps="results/hpsmc/hpsmc-plot/{pop1}---x---{pop2}.eps",
        par="results/hpsmc/hpsmc-plot/{pop1}---x---{pop2}.par"
    conda:
        "../envs/hpsmc.yaml"
    log:
        "results/logs/hpsmc/hpsmc-plot/{pop1}---x---{pop2}.log"
    benchmark:
        "results/benchmarks/hpsmc/hpsmc-plot/{pop1}---x---{pop2}.bmk"
    shell:
        "psmc_plot.pl -u 8.0e-09 -g 3 {output.eps} {input} 2> {log}"

## 4. Run simulations of divergence without post-divergence migration to compare to the original hPSMC plot
#rule simulate_hpsmc_divergence:
#    input:
#    output:
#    conda:
#    log:
#    benchmark:
#    shell:
#        " python workflow/scripts/hPSMC/hPSMC_quantify_split_time.py --Ne=? -l ? -u ? -s "
 	
    #OPTIONS:
 	#-h print help message with user options
 	#-o OUT, --out=OUT     output directory for simulations and prefix all files for the run, default="./hPSMC_sim_"
 	#-N NE, --Ne=NE        The ancestral population size to simulate, default=10,000
 	#-l LOWER, --lower=LOWER		lower bound for simulations, the most recent divergence time to be simulated
 	#-u UPPER, --upper=UPPER		upper bound for simulations, the most ancient divergence time to be simulated
 	#-s SIMS, --sims=SIMS  the number of simulations to conduct, simulations will evenly split between high and low, minimum value=2, minimum meaningful value=3
 	#-p PAR, --parallel=PAR		Number of simulations to run simultaneously
 	#-P PSMC, --PSMC=PSMC  If the psmc executable is not in your path give it's location, default = "psmc"
 	#-m MS, --ms=MS        If the ms executable is not in your path give it's location, default = "ms"
 	#-H HPSMC, --hPSMC=HPSMC
 	#If the hPSMC directory is not in your path give it's location, NOTE:  Just the directory not the script.  default = "./"

## Plot simulations using either method in step 3, with the orignal data to show the divergence between samples. 
# A) Compare simulations' pre-divergence Ne to your data, if they converge the simulations are appropriate, 
# if not reestimate Ne and repeat steps 3-5. 
# B) Identify the range of values for divergence time that intersect your hPSMC plot between 1.5 and 10 times 
# pre-divergence Ne. These are your simulations consistent with data. Your diverence time estimate is the 
# narrowest range of inconsistent simulations surrounding the consistent simulations.