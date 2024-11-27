## This rule runs model finder only (compared to running model finder and immediately going into the tree construction)
# I optionally constrained model finder (-mset) to test only those models with ascertainment bias correction (+ASC) 
# bc our data is only SNPs, and -st tells IQtree the sequence type is DNA bc it said it detected protein sequences 
rule find_model_iqtree:
    input:
        "results/plink/phylip/aut-snps-{maf}.phy",
    output:
        prefix="results/tree/aut-snps-{maf}-model-finder",
    conda:
        "../envs/iqtree2.yaml"
    log:
        "results/logs/tree/aut-snps-{maf}-model-finder.log",
    benchmark:
        "results/benchmarks/tree/aut-snps-{maf}-model-finder.bmk",
    shell:
        " iqtree2 -s {input} -st DNA -m MF -mset GTR+ASC "
        " --prefix {output.prefix} 2> {log} "

## This rule runs IQTree2 on the Phylip file from the filtered BCF created in plink.smk 
# by default, IQTree2 runs -m MFP which runs model finder plus and uses the model with the smallest BIC score to make an ML tree
# the -bb 1000 flag runs UFBoot for 1000 bootstrap replicates to generate branch support files
# the -alrt 1000 flag runs the SH-like approximate likelihood ratio test (Guindon et al., 2010) 
# for 1000 boostrap replicates which is the recommended minimum 
# both flags (-alrt, -B) will give both SH-aLRT and UFBoot support values for each branch
# when running a model with +ASC, an error will be thrown if any invariable sites are left in the alignment but IQTree(v>1.5)
# is kind enough to spit out a .varsite.phy file in this case that can be used to rerun  
rule make_iqtree_MAC1:
    input:
        "results/tree/MAC1/aut-bisnps-no5indel-nooutlier-MAC1-tree.varsites.phy",
    output:
        prefix="results/tree/MAC1/aut-bisnps-no5indel-nooutlier-MAC1-tree",
    conda:
        "../envs/iqtree2.yaml"
    log:
        "results/logs/tree/MAC1/aut-bisnps-no5indel-nooutlier-MAC1-tree.log",
    benchmark:
        "results/benchmarks/tree/MAC1/aut-bisnps-no5indel-nooutlier-MAC1-tree.bmk",
    resources:
        mem_mb=192000,
        cpus=4,
        time="23:59:59"
    shell:
        " iqtree2 -s {input} -st DNA -bb 1000 -m GTR+I+G+ASC -nt AUTO "
        " --prefix {output.prefix} 2> {log} "

rule make_iqtree_MAC3:
    input:
        "results/tree/MAC3/aut-bisnps-no5indel-nooutlier-MAC3-tree.varsites.phy",
    output:
        prefix="results/tree/MAC3/aut-bisnps-no5indel-nooutlier-MAC3-tree",
    conda:
        "../envs/iqtree2.yaml"
    log:
        "results/logs/tree/MAC3/aut-bisnps-no5indel-nooutlier-MAC3-tree.log",
    benchmark:
        "results/benchmarks/tree/MAC3/aut-bisnps-no5indel-nooutlier-MAC3-tree.bmk",
    resources:
        mem_mb=192000,
        cpus=4,
        time="23:59:59"
    shell:
        " iqtree2 -s {input} -st DNA -bb 1000 -m GTR+I+G+ASC -nt AUTO "
        " --prefix {output.prefix} 2> {log} "

rule make_iqtree_MAC5:
    input:
        "results/tree/MAC5/aut-bisnps-no5indel-nooutlier-MAC5-tree.varsites.phy",
    output:
        prefix="results/tree/MAC5/aut-bisnps-no5indel-nooutlier-MAC5-tree",
    conda:
        "../envs/iqtree2.yaml"
    log:
        "results/logs/tree/MAC5/aut-bisnps-no5indel-nooutlier-MAC5-tree.log",
    benchmark:
        "results/benchmarks/tree/MAC5/aut-bisnps-no5indel-nooutlier-MAC5-tree.bmk",
    resources:
        mem_mb=192000,
        cpus=4,
        time="23:59:59"
    shell:
        " iqtree2 -s {input} -st DNA -bb 1000 -m GTR+I+G+ASC -nt AUTO "
        " --prefix {output.prefix} 2> {log} "


## no MAC filter
#rule make_iqtree:
#    input:
#        "results/tree/aut-bisnps-no5indel-tree.varsites.phy",
#    output:
#        prefix="results/tree/aut-bisnps-no5indel-tree",
#    conda:
#        "../envs/iqtree2.yaml"
#    log:
#        "results/logs/tree/aut-bisnps-no5indel-tree.log",
#    benchmark:
#        "results/benchmarks/tree/aut-bisnps-no5indel-tree.bmk",
#    resources:
#        mem_mb=192000,
#        cpus=4,
#        time="23:59:59"
#    shell:
#        " iqtree2 -s {input} -st DNA -bb 1000 -m GTR+I+G+ASC -nt AUTO "
#        " --prefix {output.prefix} 2> {log} "