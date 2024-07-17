## this is a rule that will build and install pcangsd into
# the active conda env. 
# then, whenever you need to use pcangsd, you call the same
# conda environment, and have the flagfile as an input
# depenency to ensure that pcangsd has already been successfully
# built into that conda env.
# the active conda env.  Yep! That works nicely. -copied from Eric's post-bcf-snakeflow
rule install_pcangsd:
    params:
        hash=config["pcangsd"]["version"],
        url=config["pcangsd"]["url"]
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
        " eval \"$(ssh-agent -s)\"; ssh-add ~/.ssh/id_ed25519 && "
        " git checkout {params.hash} && "
        " pip3 install .  ) > {log} 2>&1  "

## The following 3 rules are copied from Eric's post-bcf workflow (bcftools_filter.smk & format.smk) with edits
# https://github.com/eriqande/mega-post-bcf-exploratory-snakeflows
# this rule takes our scatter intervals file and creates a region file needed for the next step
# it uses a modified version of Eric's get_scaff_members.awk script
rule make_scat_regions:
    input:
        scat_path="results/scatter_config/scatters_1200000.tsv"
    params:
        scat="{scatter}"
    output:
        "results/pca/scat_regions/{scatter}.scat_regions.tsv",
    log:
        "results/logs/pca/scat_regions/{scatter}.scat_regions.log"
    benchmark:
        "results/benchmarks/pca/scat_regions/{scatter}.scat_regions.bmk"
    shell:
        " (awk -v scat='{params.scat}' -f workflow/scripts/pca/get_scat_regions.awk {input.scat_path} > {output}) 2> {log} "


# makes beagle GL file from the PL field in BCF file.  Note that we capitate them
# and keep them around in case the individual sections are useful down the road
rule bcf2beagle_gl_scatter:
    input:
        bcf="results/bcf/all.bcf",
        csi="results/bcf/all.bcf.csi",
        regions="results/scat_regions/{scatter}.scat_regions.tsv",
        sfile=get_samples_txt
    output:
        body=temp("results/bcf/beagle-gl/sections/{scatter}.body.gz"),
        top_row=temp("results/bcf/beagle-gl/sections/{scatter}.toprow.gz"),
        beag="results/bcf/beagle-gl/sections/{scatter}.beagle-gl.gz"
    log:
        "results/logs/bcf2beagle_gl_scatter/bcf/{scatter}.log"
    benchmark:
        "results/benchmarks/bcf2beagle_gl_scatter/bcf/{scatter}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        " ( " 
        " awk -f workflow/scripts/pca/beagle3header.awk {input.sfile} | gzip -c > {output.top_row}  && "
        " bcftools view -Ou -R {input.regions} {input.bcf} |  "
        " bcftools query -f '%CHROM:%POS\\t%REF\\t%ALT[\\t%PL]\\n' | "
        " awk -f workflow/scripts/pca/pl2gl.awk | gzip -c  >  {output.body} && "
        " cat {output.top_row} {output.body} > {output.beag}  " 
        " ) 2> {log}  "



rule bcf2beagle_gl_gather:
    input: 
        header=expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_{{thin_int}}_{{thin_start}}/beagle-gl/sections/{sg}.toprow.gz", sg=first_scatter_id),
        scaff_gzs=expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_{{thin_int}}_{{thin_start}}/beagle-gl/sections/{sg}.body.gz", sg=unique_scatters)
    output:
        "results/bcf/beagle-gl/beagle-gl.gz"
    log:
        "results/logs/bcf2beagle_gl_gather/bcf/main.log"
    benchmark:
        "results/logs/bcf2beagle_gl_gather/bcf/main.bmk"
    shell:
        "cat {input.header} {input.scaff_gzs} > {output} 2> {log} "




### These are rules to get PCA from BAM files
## this rule takes the bam.filelist and generates a beagle file for use in pcangsd
#rule get_beagle_input:
#    input:
#        expand("results/angsd_bams/overlap_clipped/{s}.bam", s=sample_list),
#    output:
#        "results/pca/get-beagle-input/all.beagle.gz"
#    conda:
#        "../envs/angsd.yaml"
#    log:
#        "results/logs/pca/get-beagle-input/all.log"
#    benchmark:
#        "results/benchmarks/pca/get-beagle-input/all.bmk"
#    shell:
#        " angsd -GL 2 -baq 2 -out {output} -nThreads 10 -doGlf 2 "
#        " -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam {input} "


## this rule runs pcangsd on the beagle input file
# with no genotype posteriors, selection, or anything else
#rule run_pcangsd:
#    input:
#        flagfile="results/flags/pcangsd_installed",
#        beagle="results/pca/get-beagle-input/all.beagle.gz"
#    output:
#        "results/pca/run-pcangsd/all.cov"
#    conda:
#        "../envs/pcangsd.yaml"
#    log:
#        "results/logs/pca/run-pcangsd/all.log"
#    benchmark:
#        "results/benchmarks/pca/run-pcangsd/all.bmk"
#    shell:
#        "pcangsd -beagle {input.beagle} -o {output} 2> {log}"