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
# into the given path and make (compile) it
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
        dir=directory("results/chromcompare"),
    conda:
        "../envs/chromcompare.yaml"
    log:
        "results/logs/install_chromcompare/log.txt"
    shell:
        "(mkdir {output.dir} && cd {output.dir} && "
        " git clone {params.url} && "
        " cd Chrom-Compare  && "
        " git checkout {params.hash} && "
        " make ) > {log} 2>&1  "


## 1. Create an hPSMC.psmcfa file for each combination of 2 samples 
# for pu2fa, -c sets the chrom name and -C 50 sets a maximum depth limit (50)
# the defaults given in the hPSMC documentation say to use -q30 and -Q60 for samtools, but that returns empty files for my dataset
rule haploidize_bam_sect:
    input:
        bam=get_hpsmc_bams_in_pop,
        ref="resources/genome/OmykA.fasta",
        dir="results/chromcompare"
    output:
        "results/hpsmc/haploidize_bam_sect/{hpsmcpops}/{hpsmcchroms}_haploidized.fa",
    conda:
        "../envs/chromcompare.yaml"
    resources:
        time="23:59:59",
    log:
        "results/logs/hpsmc/haploidize-bam-sect/{hpsmcpops}/{hpsmcchroms}.log",
    benchmark:
        "results/benchmarks/hpsmc/haploidize-bam-sect/{hpsmcpops}/{hpsmcchroms}.bmk",
    shell:
        " samtools mpileup -s -f {input.ref} -q 30 -Q 30 -r {wildcards.hpsmcchroms} {input.bam} | "
        " {input.dir}/Chrom-Compare/pu2fa -C 50 -c {wildcards.hpsmcchroms} > {output} 2> {log} "


rule concat_haploidized_bam:
    input:
        expand("results/hpsmc/haploidize_bam_sect/{{hpsmcpops}}/{c}_haploidized.fa", c=hpsmcchroms),
    output:
        "results/hpsmc/haploidized_bam/{hpsmcpops}_haploidized.fa",
    log:
        "results/logs/hpsmc/concat_haploidized_bam/{hpsmcpops}.log",
    benchmark:
        "results/benchmarks/hpsmc/concat_haploidized_bam/{hpsmcpops}.log",
    shell:
        " cat {input} > {output} 2> {log} "


rule psmcfa_from_2_fastas:
    input:
        pop1="results/hpsmc/haploidized_bam/{pop1}_haploidized.fa",
        pop2="results/hpsmc/haploidized_bam/{pop2}_haploidized.fa"
    output:
        "results/hpsmc/psmcfa-from-2-fastas/{pop1}---x---{pop2}.psmcfa"
    conda:
        "../envs/hpsmc.yaml"
    resources:
        time="12:00:00",
        mem_mb=9400,
        cpus=2,
    log:
        "results/logs/hpsmc/psmcfa-from-2-fastas/{pop1}---x---{pop2}.log"
    benchmark:
        "results/benchmarks/hpsmc/psmcfa-from-2-fastas/{pop1}---x---{pop2}.bmk"
    shell:
        "python workflow/scripts/hPSMC/psmcfa_from_2_fastas.py -b100 -m50 {input.pop1} {input.pop2} > {output} 2> {log}"


## 2. run PSMC on each of the pop1---x---pop2.psmcfa files
# using same defaults as psmc
rule run_hpsmc:
    input:
        "results/hpsmc/psmcfa-from-2-fastas/{pop1}---x---{pop2}.psmcfa"
    output:
        "results/hpsmc/run-hpsmc/{pop1}---x---{pop2}.psmc"
    conda:
        "../envs/hpsmc.yaml"
    resources:
        time="23:59:59",
        mem_mb=112200,
    log:
        "results/logs/hpsmc/run-hpsmc/{pop1}---x---{pop2}.log"
    benchmark:
        "results/benchmarks/hpsmc/run-hpsmc/{pop1}---x---{pop2}.bmk"
    shell:
        "psmc -N25 -t10 -r5 -p '10+6*2+18*1+8*2+8*1' -o {output} {input} 2> {log}"

## 3. visualize hPSMC plots (using PSMC) and esimate pre-divergence Ne & upper and lower divergence time by looking at plots
## rule to plot hpsmc to visualize result
# -u [per-generation mutation rate] from https://doi.org/10.1371/journal.pgen.1010918
# -g [generation time in years] 
# -Y [maximum pop size] (in x10^4)
rule hpsmc_plot:
    input:
        "results/hpsmc/run-hpsmc/{pop1}---x---{pop2}.psmc"
    output:
        "results/hpsmc/hpsmc-plot/{pop1}---x---{pop2}",
        #par="results/hpsmc/hpsmc-plot/{pop1}---x---{pop2}.par"
    conda:
        "../envs/hpsmc.yaml"
    log:
        "results/logs/hpsmc/hpsmc-plot/{pop1}---x---{pop2}.log"
    benchmark:
        "results/benchmarks/hpsmc/hpsmc-plot/{pop1}---x---{pop2}.bmk"
    shell:
        "psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" -Y 40 {output} {input} 2> {log}"


rule hpsmc_plot_multiple:
    input:
        "results/hpsmc/run-hpsmc/w_fk_boulder---x---nanita.psmc",
        "results/hpsmc/run-hpsmc/w_fk_boulder---x---steelman.psmc",
        "results/hpsmc/run-hpsmc/w_fk_boulder---x---san_juan.psmc",
        "results/hpsmc/run-hpsmc/w_fk_boulder---x---rio_grande.psmc",
        "results/hpsmc/run-hpsmc/w_fk_boulder---x---como.psmc",
        "results/hpsmc/run-hpsmc/w_fk_boulder---x---s_hayden.psmc",
        "results/hpsmc/run-hpsmc/greenback---x---w_fk_boulder.psmc",
    #params:
    #    get_comma_sep_hpsmc_names
    output:
        "results/hpsmc/hpsmc-plot/w_fk_boulder-x-all",
        #par="results/hpsmc/hpsmc-plot/w_fk_boulder-x-all.par"
    conda:
        "../envs/hpsmc.yaml"
    log:
        "results/logs/hpsmc/hpsmc-plot/w_fk_boulder-x-all.log"
    benchmark:
        "results/benchmarks/hpsmc/hpsmc-plot/w_fk_boulder-x-all.bmk"
    shell:
        "psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" "
        "-M w_fk_boulder---x---nanita,w_fk_boulder---x---steelman,w_fk_boulder---x---san_juan,w_fk_boulder---x---rio_grande,w_fk_boulder---x---como,w_fk_boulder---x---s_hayden,greenback---x---w_fk_boulder"
        "{output} {input} 2> {log}"



## 4. Run simulations of divergence without post-divergence migration to compare to the original hPSMC plot
rule get_hpsmc_divergence_sh:
    input:
        "results/hpsmc/run-hpsmc/greenback---x---s_hayden.psmc",
    output:
        sh="results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_quantify_split_time.sh"
    params:
        pfx="results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_",
    conda:
        "../envs/hpsmc-split.yaml"
    log:
        "results/logs/hpsmc/split-time-sim/greenback---x---s_hayden_split_sim.log"
    benchmark:
        "results/benchmarks/hpsmc/split-time-sim/greenback---x---s_hayden_split_sim.bmk"
    shell:
        " python workflow/scripts/hPSMC/quantify_split_time.py --Ne=30000 -l 0 -u 350000 -s 10 -p 10 "
        " -o {params.pfx} --hPSMC workflow/scripts/hPSMC/ {input} > {output.sh} 2> {log} "

rule simulate_hpsmc_divergence:
    input:
        sh="results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_quantify_split_time.sh",
    output:
        log="results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_split_result.log"
    resources:
        time="23:59:59",
        mem_mb=112200,
    conda:
        "../envs/hpsmc-split.yaml"
    benchmark:
        "results/benchmarks/hpsmc/split-time-sim/greenback---x---s_hayden_split_result.bmk"
    shell:
        " chmod +rwx {input.sh} && "
        " bash ./{input.sh} 2> {output.log} "
 	
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

rule hpsmc_plot_sims:
    input:
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_0.ms_sim.psmc",
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_38888.ms_sim.psmc",
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_77777.ms_sim.psmc",
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_116666.ms_sim.psmc",
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_155555.ms_sim.psmc",
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_194444.ms_sim.psmc",
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_233333.ms_sim.psmc",
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_272222.ms_sim.psmc",
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_311111.ms_sim.psmc",
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden_350000.ms_sim.psmc",
    #params:
    #    get_comma_sep_hpsmc_names
    output:
        "results/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden",
        #par="results/hpsmc/hpsmc-plot/w_fk_boulder-x-all.par"
    conda:
        "../envs/hpsmc.yaml"
    log:
        "results/logs/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden-plot.log"
    benchmark:
        "results/benchmarks/hpsmc/split-time-sim/greenback---x---s_hayden/greenback---x---s_hayden-plot.bmk"
    shell:
        " psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" "
        " -M 0,38888,77777,116666,155555,194444,233333,272222,311111,350000 "
        " {output} {input} 2> {log} "