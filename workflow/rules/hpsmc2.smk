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
        "../envs/chromcompare.yaml"
    log:
        "results/logs/hpsmc-test/haploidize-bam-sect/{hpsmcpops}/{hpsmcchroms}.log",
    benchmark:
        "results/benchmarks/hpsmc-test/haploidize-bam-sect/{hpsmcpops}/{hpsmcchroms}.bmk",
    shell:
        " samtools mpileup -s -f {input.ref} -q 30 -Q 30 -r {wildcards.hpsmcchroms} {input.bam} | "
        " {input.dir}/Chrom-Compare/pu2fa -C 50 -c {wildcards.hpsmcchroms} > {output} 2> {log} "


rule test_concat_haploidized_bam:
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


rule test_psmcfa_from_2_fastas:
    input:
        pop1="results/hpsmc-test/haploidized_bam/{pop1}_haploidized.fa",
        pop2="results/hpsmc-test/haploidized_bam/{pop2}_haploidized.fa"
    output:
        "results/hpsmc-test/psmcfa-from-2-fastas/{pop1}---x---{pop2}.psmcfa"
    conda:
        "../envs/hpsmc.yaml"
    resources:
        time="12:00:00",
        mem_mb=9400,
        cpus=2,
    log:
        "results/logs/hpsmc-test/psmcfa-from-2-fastas/{pop1}---x---{pop2}.log"
    benchmark:
        "results/benchmarks/hpsmc-test/psmcfa-from-2-fastas/{pop1}---x---{pop2}.bmk"
    shell:
        "python workflow/scripts/hPSMC/JO_psmcfa_from_2_fastas.py -b5 -m2 {input.pop1} {input.pop2} > {output} 2> {log}"

## 2. run PSMC on each of the pop1---x---pop2.psmcfa files
# using same defaults as psmc
rule test_run_hpsmc:
    input:
        "results/hpsmc-test/psmcfa-from-2-fastas/{pop1}---x---{pop2}.psmcfa"
    output:
        "results/hpsmc-test/run-hpsmc/{pop1}---x---{pop2}.psmc"
    conda:
        "../envs/hpsmc.yaml"
    resources:
        time="23:59:59",
        mem_mb=112200,
    log:
        "results/logs/hpsmc-test/run-hpsmc/{pop1}---x---{pop2}.log"
    benchmark:
        "results/benchmarks/hpsmc-test/run-hpsmc/{pop1}---x---{pop2}.bmk"
    shell:
        "psmc -N25 -t10 -r5 -p '10+6*2+18*1+8*2+8*1' -o {output} {input} 2> {log}"

## 3. visualize hPSMC plots (using PSMC) and esimate pre-divergence Ne & upper and lower divergence time by looking at plots
## rule to plot hpsmc to visualize result
# -u [per-generation mutation rate] from https://doi.org/10.1371/journal.pgen.1010918
# -g [generation time in years] 
rule test_hpsmc_plot:
    input:
        "results/hpsmc-test/run-hpsmc/{pop1}---x---{pop2}.psmc"
    output:
        "results/hpsmc-test/hpsmc-plot/{pop1}---x---{pop2}",
        #par="results/hpsmc/hpsmc-plot/{pop1}---x---{pop2}.par"
    conda:
        "../envs/hpsmc.yaml"
    log:
        "results/logs/hpsmc-test/hpsmc-plot/{pop1}---x---{pop2}.log"
    benchmark:
        "results/benchmarks/hpsmc-test/hpsmc-plot/{pop1}---x---{pop2}.bmk"
    shell:
        "psmc_plot.pl -u 8.0e-09 -g 3 -P \"below\" {output} {input} 2> {log}"