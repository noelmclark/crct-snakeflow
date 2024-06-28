# copied from Eric's calling.smk

## the following 3 rules are used to create the scaffold group, chromo, and scatter interval lists
# which are used to joint call sample by chromo or scaffold group using the genomics db 
rule make_scaff_group_interval_lists:
    input:
        scaff_groups = config["scaffold_groups"]
    output:
        "results/calling/interval_lists/{scaff_group}.list"
    log:
        "results/logs/calling/make_interval_lists/{scaff_group}.log"
    shell:
        " awk -v sg={wildcards.scaff_group} 'NR>1 && $1 == sg {{print $2}}' {input.scaff_groups} > {output} 2> {log};"

rule make_chromo_interval_lists:
    output:
        "results/calling/interval_lists/{chromo}.list"
    log:
        "results/logs/calling/make_interval_lists/{chromo}.log"
    shell:
        " echo {wildcards.chromo} > {output} 2> {log};"

rule make_scatter_interval_lists:
    input:
        scatters_file= config["scatter_intervals_file"]
    log:
        "results/logs/calling/make_scatter_interval_lists/{sg_or_chrom}/{scatter}.log"
    output:
        "results/calling/scatter_interval_lists/{sg_or_chrom}/{scatter}.list"
    shell:
        " awk -v sgc={wildcards.sg_or_chrom} -v scat={wildcards.scatter} ' "
        "    NR>1 && $1 == sgc && $2==scat {{printf(\"%s:%s-%s\\n\", $3, $4, $5)}} "
        " ' {input.scatters_file} > {output} 2> {log};"


## the following 2 rules create gvcf files for each sample to be used in the genomics db
# the first rule creates these per sample by sections of scaffold group or chromo
# the second rule concats these sections into one gvcf per sample
rule make_gvcf_sections:
    input:
        bam="results/mapping/gatk-rmdup/{sample}.bam",
        bai="results/mapping/gatk-rmdup/{sample}.bai",
        ref="resources/genome/OmykA.fasta",
        idx="resources/genome/OmykA.dict",
        fai="resources/genome/OmykA.fasta.fai",
        interval_list="results/calling/interval_lists/{sg_or_chrom}.list"
    output:
        gvcf="results/calling/gvcf_sections/{sample}/{sg_or_chrom}.g.vcf.gz",
        idx="results/calling/gvcf_sections/{sample}/{sg_or_chrom}.g.vcf.gz.tbi",
    conda:
        "../envs/gatk.yaml"
    log:
        stderr="results/logs/calling/make_gvcf_sections/{sample}/{sg_or_chrom}.stderr",
        stdout="results/logs/calling/make_gvcf_sections/{sample}/{sg_or_chrom}.stdout",
    benchmark:
        "results/benchmarks/calling/make_gvcfs_sections/{sample}/{sg_or_chrom}.bmk"
    params:
        java_opts="-Xmx4g",
        conf_pars=config["params"]["gatk"]["HaplotypeCaller"]
    resources:
        time="1-00:00:00",
        mem_mb = 4600,
        cpus = 1
    threads: 1
    shell:
        "gatk --java-options \"{params.java_opts}\" HaplotypeCaller "
        " -R {input.ref} "
        " -I {input.bam} "
        " -O {output.gvcf} "
        " -L {input.interval_list} "
        " --native-pair-hmm-threads {threads} "
        " {params.conf_pars} "
        " -ERC GVCF > {log.stdout} 2> {log.stderr} "


## This makes a single GVCF file per individual sample. 
rule concat_gvcf_sections:
    input: 
        expand("results/calling/gvcf_sections/{s}/{sgc}.g.vcf.gz", s=sample_list, sgc=sg_or_chrom)
    output:
        gvcf="results/calling/gvcf/{sample}.g.vcf.gz",
        idx="results/calling/gvcf/{sample}.g.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/calling/concat_gvcf_sections/{sample}.txt"
    benchmark:
        "results/benchmarks/calling/concat_gvcf_sections/{sample}.bmk"
    params:
        opts=" --naive "
    shell:
        " bcftools concat {params.opts} -O z {input} > {output.gvcf} 2>{log}; "
        " bcftools index -t {output.gvcf} "



## The following 2 rules create the genomics db by chromo and scaff_groups 
# which are used to joint call variants for all sample by section specified
# From https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport 
# GenomicsDBImport uses temporary disk storage during import. The amount of temporary disk storage required 
# can exceed the space available, especially when specifying a large number of intervals. 
# The command line argument `--tmp-dir` can be used to specify an alternate temporary storage location with sufficient space.

rule import_genomics_db_by_chromo:
    input:
        gvcfs=expand("results/calling/gvcf_sections/{s}/{chromo}.g.vcf.gz", s=sample_list, chromo=unique_chromosomes),
        gvcf_idx=expand("results/calling/gvcf_sections/{s}/{chromo}.g.vcf.gz.tbi", s=sample_list, chromo=unique_chromosomes),
    output:
        gdb=directory("results/calling/genomics_db/{chromo}")
    conda:
        "../envs/gatk.yaml"
    log:
        "results/logs/calling/import_genomics_db/{chromo}.log"
    benchmark:
        "results/benchmarks/calling/import_genomics_db/{chromo}.bmk"
    params:
        java_opts="-Xmx4g -Xms4g",
        my_opts=chromo_import_gdb_opts,
        extras="--tmp-dir=results/snake-tmp",
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    shell:
        " VS=$(for i in {input.gvcfs}; do echo -V $i; done); "  # make a string like -V file1 -V file2
        " gatk --java-options {params.java_opts} GenomicsDBImport "
        " $VS "
        " {params.my_opts} {params.extras}"
        " --genomicsdb-workspace-path {output.gdb} "
        " 2> {log} "


rule import_genomics_db_by_scaffold_group:
    input:
        gvcfs=expand("results/calling/gvcf_sections/{s}/{scaff_group}.g.vcf.gz", s=sample_list, scaff_group=unique_scaff_groups),
        gvcf_idx=expand("results/calling/gvcf_sections/{s}/{scaff_group}.g.vcf.gz.tbi", s=sample_list, scaff_group=unique_scaff_groups),
    output:
        gdb=directory("results/calling/genomics_db/{scaff_group}"),
    conda:
        "../envs/gatk.yaml"
    log:
        "results/logs/calling/import_genomics_db/{scaff_group}.log"
    benchmark:
        "results/benchmarks/calling/import_genomics_db/{scaff_group}.bmk"
    params:
        java_opts="-Xmx4g -Xms4g",
        my_opts=scaff_group_import_gdb_opts,
        extras="--tmp-dir=results/snake-tmp",
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    shell:
        " VS=$(for i in {input.gvcfs}; do echo -V $i; done); "  # make a string like -V file1 -V file2
        " gatk --java-options {params.java_opts} GenomicsDBImport "
        " $VS "
        " {params.my_opts} {params.extras}"
        " --genomicsdb-workspace-path {output.gdb} "
        " 2> {log} "


## The next rule uses GenotypeGVCFs to do joint genotyping using a genomics db 
# and a list of smaller pieces of the chroms and scaffold groups (scatters) 
# to get one vcf file per chrom or scaff group with all of the samples in it

rule vcf_scattered_from_gdb:
    input:
        gdb="results/calling/genomics_db/{sg_or_chrom}",
        scatters="results/calling/scatter_interval_lists/{sg_or_chrom}/{scatter}.list",
        ref="resources/genome/OmykA.fasta",
        fai="resources/genome/OmykA.fasta.fai",
        idx="resources/genome/OmykA.dict",
    output:
        vcf="results/calling/vcf_sections/{sg_or_chrom}/{scatter}.vcf.gz",
        idx="results/calling/vcf_sections/{sg_or_chrom}/{scatter}.vcf.gz.tbi",
    conda:
        "../envs/gatk.yaml"
    log:
        "results/logs/calling/vcf_scattered_from_gdb/{sg_or_chrom}/{scatter}.txt"
    benchmark:
        "results/benchmarks/calling/vcf_scattered_from_gdb/{sg_or_chrom}/{scatter}.bmk"
    params:
        java_opts="-Xmx4g",
        extra=" --genomicsdb-shared-posixfs-optimizations --only-output-calls-starting-in-intervals " #from Eric, idk meaning
    resources:
        mem_mb = 11750,
        cpus = 2,
        time = "23:59:59"
    threads: 2
    shell:
        " gatk --java-options {params.java_opts} GenotypeGVCFs "
        "  {params.extra} "
        "  -L {input.scatters} "
        "  -R {input.ref}  "
        "  -V gendb://{input.gdb} "
        "  -O {output.vcf} 2> {log} "


## This rule takes the vcf files for each small chunk (scatter) of the chroms and scaffold groups
# and concats them back together based on chrom or scaffold group. So we end up with one vcf file
# per chrom or scaffold group that contains all the samples variant info. 
rule gather_scattered_vcfs:
    input:
        vcf=lambda wc: get_scattered_vcfs(wc, ""),
        tbi=lambda wc: get_scattered_vcfs(wc, ".tbi"),
    output:
        vcf="results/calling/vcf_sections/{sg_or_chrom}.vcf.gz",
        tbi="results/calling/vcf_sections/{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/calling/gather_scattered_vcfs/{sg_or_chrom}.txt"
    benchmark:
        "results/benchmarks/calling/gather_scattered_vcfs/{sg_or_chrom}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Oz {input.vcf} > {output.vcf}; "
        " bcftools index -t {output.vcf})  2>{log}; "


## From Eric's workflow
# this is a little rule we throw in here so that we can mark
# an individual as missing data (./. or .|.) when it has a read
# depth of 0, because GATK now marks those as 0/0,
# see https://gatk.broadinstitute.org/hc/en-us/community/posts/4476803114779-GenotypeGVCFs-Output-no-call-as-reference-genotypes?page=1#community_comment_6006727219867
# this also adds an INFO field NMISS, which gives the number of samples missing a call.
# 8/26/22: This has been updated to also mark genotypes as missing if they have a PL of 0,0,0.
rule correct_missing_vcf_sect:
    input:
        vcf="results/calling/vcf_sections/{sg_or_chrom}.vcf.gz"
    output:
        vcf="results/calling/corrected_missing_vcf_sect/{sg_or_chrom}.vcf.gz",
        tbi="results/calling/corrected_missing_vcf_sect/{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/calling/correct_missing_vcf_sect/{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/calling/correct_missing_vcf_sect/{sg_or_chrom}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools +setGT {input.vcf} -- -t q -n . -i 'FMT/DP=0 | (FMT/PL[:0]=0 & FMT/PL[:1]=0 & FMT/PL[:2]=0)' | "
        " bcftools +fill-tags - -- -t 'NMISS=N_MISSING' | "
        " bcftools view -Oz - > {output.vcf}; "
        " bcftools index -t {output.vcf}) 2> {log} "



