## this is a rule that will build and install pcangsd into
# the active conda env. 
# then, whenever you need to use pcangsd, you call the same
# conda environment, and have the flagfile as an input
# depenency to ensure that pcangsd has already been successfully
# built into that conda env.
# the active conda env.  Yep! That works nicely. -copied from Eric's post-bcf-snakeflow
rule install_pcangsd:
    params:
        url=config["params"]["pcangsd"]["url"]
    output:  
        flagfile=touch("results/flags/pcangsd_installed")
    conda:
        "../envs/pcangsd.yaml"
    log:
        "results/logs/install_pcangsd/log.txt"
    shell:
        "(TMP=$(mktemp -d) && cd $TMP && "
        " git clone {params.url} && "
        " cd pcangsd  && "
        " python setup.py build_ext --inplace && "
        " pip3 install -e .  ) > {log} 2>&1  "


## this rule takes the bam.filelist and generates a beagle file for use in pcangsd
# the {post_id} wildcard determines which bam files are used as input - ALL is all samples,
# CRCT is just the blue/green lineage individuals & OUT is just the outgroups
rule get_beagle_input:
    input:
        get_pcangsd_bam_filelist
    output:
        "results/pca/get-beagle-input/{post_id}.beagle.gz"
    conda:
        "../envs/angsd.yaml"
    log:
        "results/logs/pca/get-beagle-input/{post_id}.log"
    benchmark:
        "results/benchmarks/pca/get-beagle-input/{post_id}.bmk"
    shell:
        " angsd -GL 2 -baq 2 -out {output} -nThreads 10 -doGlf 2 "
        " -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam {input} "


## this rule runs pcangsd on the beagle input file
# with no genotype posteriors, selection, or anything else
rule run_pcangsd:
    input:
        flagfile="results/flags/pcangsd_installed",
        beagle="results/pca/get-beagle-input/{post_id}.beagle.gz"
    output:
        "results/pca/run-pcangsd/{post_id}.cov"
    conda:
        "../envs/pcangsd.yaml"
    log:
        "results/logs/pca/run-pcangsd/{post_id}.log"
    benchmark:
        "results/benchmarks/pca/run-pcangsd/{post_id}.bmk"
    shell:
        "pcangsd -beagle {input.beagle} -o {output} 2> {log}"