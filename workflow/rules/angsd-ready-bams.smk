## These are rules that take the mapped mkdup BAMs, and the results
# of running GATK to get the indels, and do:
#    1. Clipping over overlaps
#    2. AddOrReplaceReadGroups 
#    3. Indel realignment
# Apparently indel realignment is a bit cumbersome and not necessary in some cases.
# An alternative presented by (https://academic.oup.com/bioinformatics/article/27/8/1157/227268) 
# is just using base alignment quality (-baq) option in ANGSD. 


## 1. Simple clipping of overlaps in the mkduped BAMs (does not remove duplicates)
rule clip_overlaps:
    input:
        "results/mapping/gatk-rmdup/{sample}.bam"
    output:
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",
        bai="results/angsd_bams/overlap_clipped/{sample}.bai"
    log:
        clip="results/logs/angsd_bams/clip_overlaps/{sample}.log",
        index="results/logs/angsd_bams/clip_overlaps/index-{sample}.log"
    conda:
        "../envs/bamutil_samtools.yaml"
    benchmark:
        "results/benchmarks/angsd_bams/clip_overlaps/{sample}.bmk"
    shell:
        " bam clipOverlap --in {input} --out {output.bam} --stats 2> {log.clip} && "
        " samtools index {output.bam} -o {output.bai} 2> {log.index}"


## 2. This is a new rule I wrote because HaplotyeCaller was throwing an error in my make_gvcf_sections rule with the bam files
# from samples that had multiple units because I erroneously had the SM field of @RG include unit. 
# They had multiple @RG tags listed so HaplotyeCaller thought there were two samples 
# in these BAM files which it doesn't like. The rule worked for the samples with only one unit though.  
rule fix_RG_sample:
    input:
        "results/mapping/mapped/{sample}---{unit}.sorted.bam",
    output:
        "results/angsd_bams/RG-fixed/{sample}---{unit}.sorted.bam",
    log:
        "results/logs/angsd_bams/RG-fixed/{sample}---{unit}.log"
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "results/benchmarks/angsd_bams/RG-fixed/{sample}---{unit}.bmk"
    params:
        RRG=replace_read_group,
        #extra="--CREATE_INDEX --TMP_DIR results/snake-tmp"
    #resources:
    #    cpus = 1,
    #    mem_mb=112200,
    #threads: 30,
    shell:
        " gatk AddOrReplaceReadGroups "
        #" {params.extra} "
        " -I {input} "
        " -O {output} "
        " {params.RRG} "
        " 2> {log} "

# this is a workaround so I can request the output from the fix_RG_sample rule in the Snakefile
rule echo_RG_fixed:
    input:
        get_RG_fixed_bams_of_common_sample
    output:
        "results/angsd_bams/echo-RG-fixed/{sample}.txt",
    log:
        "results/logs/angsd_bams/echo-RG-fixed/{sample}.log"
    shell:
        " echo {input} > {output} 2> {log} "


## 3. Indel realignment - I don't actually use this in my workflow - I use the angsd -baq 2 option instead per Eric's counsel
# See Eric's full snakeflow for these rules