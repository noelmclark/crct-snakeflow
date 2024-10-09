## This rule runs IQTree2 on the Phylip file from the filtered BCF created in plink.smk 
# by default, IQTree2 runs -m MFP which runs model finder plus 
# and uses the model with the smallest BIC score to make an ML tree
# the -B 1000 flag runs UFBoot for 1000 bootstrap replicates to generate branch support files
# the -alrt 1000 flag runs the SH-like approximate likelihood ratio test (Guindon et al., 2010) 
# for 1000 boostrap replicates which is the recommended minimum 
# both flags (-alrt, -B) will give both SH-aLRT and UFBoot support values for each branch  
rule make_iqtree:
    input:
        "results/plink/phylip/all-samp-no-y.phy",
    output:
        prefix="results/tree/all-samp-no-y",
    conda:
        "../envs/iqtree2.yaml"
    log:
        "results/logs/tree/all-samp-no-y.log",
    benchmark:
        "results/benchmarks/tree/all-samp-no-y.bmk",
    resources:
        mem_mb=11750,
        cpus=2,
        time="23:59:59"
    shell:
        " iqtree2 -s {input} -alrt 1000 -B 1000 -T AUTO "
        " --prefix {output.prefix} 2> {log} "