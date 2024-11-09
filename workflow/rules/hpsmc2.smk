

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