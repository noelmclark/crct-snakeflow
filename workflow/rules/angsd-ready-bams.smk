## These are rules that take the mapped mkdup BAMs, and the results
# of running GATK to get the indels, and do:
#    1. Clipping over overlaps
#    2. Indel realignment
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


## Rule to validate BAM files to see if that is the reason why some of them are failing out of the RG fix rule...
rule validate_bam:
    input:
        bam="results/mapping/gatk-rmdup/{sample}.bam",
    output:
        "results/angsd_bams/rmdup-validated/{sample}.bam",
    log:
        "results/logs/angsd_bams/validated/rmdup-{sample}.log"
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "results/benchmarks/angsd_bams/validated/rmdup-{sample}.bmk"
    shell:
        " gatk ValidateSamFile "
        " -I {input.bam} "
        " --MODE SUMMARY "
        " > {output} 2> {log} " 


## This is a new rule I wrote because HaplotyeCaller was throwing an error in my make_gvcf_sections rule with the bam files
# from samples that had multiple units. They had multiple @RG tags listed so HaplotyeCaller thought there were two samples 
# in these BAM files which it doesn't like. The rule worked for the samples with only one unit.  
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

rule echo_RG_fixed:
    input:
        get_RG_fixed_bams_of_common_sample
    output:
        "results/angsd_bams/echo-RG-fixed/{sample}.txt",
    log:
        "results/logs/angsd_bams/echo-RG-fixed/{sample}.log"
    shell:
        " echo {input} > {output} 2> {log} "

#rule index_RG_sample:
#    input:
#        bam="results/angsd_bams/RG-fixed/{sample}.bam"
#    output:
#        bai="results/angsd_bams/RG-fixed/{sample}.bai"
#    log:
#        "results/logs/angsd_bams/RG-fixed/index-{sample}.log"
#    conda:
#        "../envs/bwa2sam.yaml"
#    benchmark:
#        "results/benchmarks/angsd_bams/RG-fixed/index-{sample}.bmk"
#    shell:
#        " samtools index {input.bam} -o {output.bai} 2> {log} "


## 2. Indel realignment
# The next chunk of rules is for doing indel realignment.
# If we want to make species specific indel VCFs we need to specify the 
# indel groups (species/lineages) we want to recognize 
rule species_sample_lists:
    params:
        sams=get_igrp_sample_names
    output:
        sample_list="results/angsd_bams/indel_realign/igrp_lists/{igrp}.args",
    log:
        "results/logs/angsd_bams/indel_realign/species_sample_lists/{igrp}.log"
    shell:
        " for i in {params.sams}; do echo $i; done > {output.sample_list} 2>{log} "



# When we have multiple species that we are dealing with,
# then we have a step of making species specific indel VCFs
rule make_species_specific_indel_vcfs:
    input:
        indel="results/hard_filtering/indels-{sg_or_chrom}.vcf.gz",
        indel_idx="results/hard_filtering/indels-{sg_or_chrom}.vcf.gz.tbi",
        sample_list="results/angsd_bams/indel_realign/igrp_lists/{igrp}.args",
        ref="resources/genome/OmykA.fasta"
    output:
         vcf=temp("results/angsd_bams/indel_realign/species-spec-indel-sections/{igrp}/{sg_or_chrom}.vcf.gz"),
         idx=temp("results/angsd_bams/indel_realign/species-spec-indel-sections/{igrp}/{sg_or_chrom}.vcf.gz.tbi"),
    log:
        "results/logs/angsd_bams/indel_realign/make_species_specific_indel_vcfs/{igrp}/{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/angsd_bams/indel_realign/make_species_specific_indel_vcfs/{igrp}/selectvariants-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    shell:
        " IGRP={wildcards.igrp};                           "  
        " if [ $IGRP = \"__ALL\" ]; then                   "   # just hard-link the files in this case
        "    (ln  {input.indel} {output.vcf};              "
        "    ln  {input.indel_idx} {output.idx}) 2>{log};  "
        " else                                        "
        "   gatk SelectVariants -R {input.ref}        "
        "     -V {input.indel}                        "
        "     -O {output.vcf}                         "
        "     --sample-name {input.sample_list}       "
        "     --exclude-non-variants 2>{log};         "
        " fi                                          "



# We need to get a VCF file holding all the known indels.
# We will use the ones that we hardfiltered.  We only keep
# the first sample in each to make it a lot smaller.  (We
# could probably make a VCF with no samples in it, but it
# is easy enough to just keep one sample, and we know that
# will work).  We do this on a chromosome or scaff group basis
rule get_known_indels:
    input:
        indel="results/angsd_bams/indel_realign/species-spec-indel-sections/{igrp}/{sg_or_chrom}.vcf.gz",
        indel_idx="results/angsd_bams/indel_realign/species-spec-indel-sections/{igrp}/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf=temp("results/angsd_bams/indel_realign/known-indel-sections/{igrp}/{sg_or_chrom}.vcf.gz"),
    log:
        "results/logs/angsd_bams/indel_realign/get_known_indels/{igrp}/{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/angsd_bams/indel_realign/get_known_indels/bcftools-{igrp}-{sg_or_chrom}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "SAMP=$(bcftools query -l {input.indel} | head -n 1); "
        " bcftools view -Oz -s $SAMP {input.indel} > {output.vcf} 2> {log}; "
        " bcftools index -t {output.vcf} 2>> {log}; "


# now we concat the above outputs all together into a single vcf for each species/lineage group
rule bcfconcat_known_indels:
    input:
        expand("results/angsd_bams/indel_realign/known-indel-sections/{{igrp}}/{sgc}.vcf.gz", sgc = sg_or_chrom)
    output:
        vcf="results/angsd_bams/indel_realign/known-indels/{igrp}/all-indels.vcf.gz",
        idx="results/angsd_bams/indel_realign/known-indels/{igrp}/all-indels.vcf.gz.tbi"
    log:
        "results/logs/angsd_bams/indel_realign/bcfconcat_known_indels/{igrp}.txt"
    benchmark:
        "results/benchmarks/angsd_bams/indel_realign/bcfconcat_known_indels/bcf_concat-{igrp}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Oz {input} > {output.vcf} 2> {log}; "
        " bcftools index -t {output.vcf})  2>> {log}; "








### The following 3 steps realign indels using GATK3 bc they got rid of IndelRealigner in GATK4...
## IDK whats going on
## actually that's a lie, I know the last step is what actually realigns the indels for each group using IndelRealigner

# Download the gatk3.8 jar and then run gatk3-register. This is necessary
# because the licensing for gatk3 doesn't allow the jar to be distributed
# directly off of conda.  Apparently this only has to be done for
# Mac OS X---gatk3-register is not even available on the Linux
# conda package.  What-evs.
rule gatk3_register:
    output:
        dir=directory("results/gatk3.8_jar"),
        jar="results/gatk3.8_jar/GenomeAnalysisTK.jar"
    log:
        "results/logs/gatk3_register/log.txt"
    benchmark:
        "results/benchmarks/gatk3_register/gatk3_register.bmk",
    conda:
        "../envs/gatk3.8.yaml"
    params:
        jopts="-Xmx4g"
    shell:
        "( wget -P {output.dir} https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2 && "
        "  bzip2 -d  {output.dir}/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2  && "
        "  tar --directory {output.dir} -xvf {output.dir}/GenomeAnalysisTK-3.8-0-ge9d806836.tar && "
        "  mv {output.dir}/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar {output.dir} && "
        "  if [ $(uname -s) = \"Darwin\" ]; then gatk3-register {output.jar}; fi"
        "  ) > {log} 2>&1 "





# Run RealignerTargetCreator from GATK 3.8.0.
rule realigner_target_creator:
    input:
        vcf="results/angsd_bams/indel_realign/known-indels/{igrp}/all-indels.vcf.gz",
        idx="results/angsd_bams/indel_realign/known-indels/{igrp}/all-indels.vcf.gz.tbi",
        jar="results/gatk3.8_jar/GenomeAnalysisTK.jar"
    output:
        "results/angsd_bams/indel_realign/realigner-targets/{igrp}/realigner-targets.intervals"   
    log:
        "results/logs/angsd_bams/indel_realign/realigner_target_creator/{igrp}.txt"
    benchmark:
        "results/benchmarks/angsd_bams/indel_realign/realigner_target_creator/{igrp}.bmk",
    conda:
        "../envs/gatk3.8.yaml"
    params:
        jopts="-Xmx4g"
    threads: 4
    shell:
        "gatk3 {params.jopts} -T RealignerTargetCreator -nt {threads} "
        " -R resources/genome.fasta "
        " -o {output} "
        " -known {input.vcf} 2> {log} "



# realign indels on each bam using the known indel variation
# stored in the realigner_target_creator step
rule indel_realigner:
    input:
        vcf="results/angsd_bams/indel_realign/known-indels/{igrp}/all-indels.vcf.gz",
        idx="results/angsd_bams/indel_realign/known-indels/{igrp}/all-indels.vcf.gz.tbi",
        jar="results/gatk3.8_jar/GenomeAnalysisTK.jar",
        intervals="results/realigner-targets/{igrp}/realigner-targets.intervals",
        fasta="resources/genome.fasta",
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",
        bai="results/angsd_bams/overlap_clipped/{sample}.bam.bai"
    output:
        bam="results/angsd_bams/indel_realigned/{igrp}/{sample}.bam"
    log:
        "results/logs/angsd_bams/indel_realign/indel_realigner/{igrp}/{sample}.txt"
    benchmark:
        "results/benchmarks/angsd_bams/indel_realign/indel_realigner/{igrp}/{sample}.bmk",
    conda:
        "../envs/gatk3.8.yaml"
    params:
        jopts="-Xmx4g"
    shell:
        "gatk3 {params.jopts}  \"-Dsamjdk.compression_level=9\" "
        " -T IndelRealigner  "
        " -R {input.fasta} "
        " -I {input.bam} "
        " -targetIntervals {input.intervals} "
        " -known {input.vcf} "
        " --consensusDeterminationModel KNOWNS_ONLY "
        " -LOD 0.4 " # this is recommended at https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-methods-and-algorithms/Local_Realignment_around_Indels.md
        " -o {output.bam} "
        " > {log} 2>&1 "