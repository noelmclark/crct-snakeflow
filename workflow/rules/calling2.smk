## something in these last 3 rules has a syntax error

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
        "results/logs/calling/vcf_scattered_from_gdb/{sg_or_chrom}.txt"
    benchmark:
        "results/benchmarks/calling/vcf_scattered_from_gdb/{sg_or_chrom}.bmk"
    params:
        java_opts="-Xmx4g"
        extra=" --genomicsdb-shared-posixfs-optimizations --only-output-calls-starting-in-intervals " #from Eric, idk meaning
    resources:
        mem_mb = 11750,
        cpus = 2,
        time = "1-00:00:00"
    threads: 2
    shell:
        " gatk --java-options {param.java_opts} GenotypeGVCFs "
        #"  {params.extra} "
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








## The next rule uses GenotypeGVCFs to do joint genotyping using a genomics db 
# and a list of smaller pieces of the chroms and scaffold groups (scatters) 
# to get one vcf file per chrom or scaff group with all of the samples in it

rule vcf_scattered_from_gdb:
    input:
        gdb="results/calling/genomics_db/{sg_or_chrom}",
        scatters=expand("results/calling/scatter_interval_lists/{sgc}/{scat}.list", sgc=sg_or_chrom, scat=unique_scats),
        ref="resources/genome/OmykA.fasta",
        fai="resources/genome/OmykA.fasta.fai",
        idx="resources/genome/OmykA.dict",
    output:
        vcf="results/calling/vcf_sections/{sg_or_chrom}/{scatter}.vcf.gz",
        idx="results/calling/vcf_sections/{sg_or_chrom}/{scatter}.vcf.gz.tbi",
    conda:
        "../envs/gatk.yaml"
    log:
        "results/logs/calling/vcf_scattered_from_gdb/{sg_or_chrom}.txt"
    benchmark:
        "results/benchmarks/calling/vcf_scattered_from_gdb/{sg_or_chrom}.bmk"
    params:
        java_opts="-Xmx4g"
        extra=" --genomicsdb-shared-posixfs-optimizations --only-output-calls-starting-in-intervals " #from Eric, idk meaning
    resources:
        mem_mb = 11750,
        cpus = 2,
        time = "1-00:00:00"
    threads: 2
    shell:
        " gatk --java-options {param.java_opts} GenotypeGVCFs "
        #"  {params.extra} "
        "  -L {input.scatters} "
        "  -R {input.ref}  "
        "  -V gendb://{input.gdb} "
        "  -O {output.vcf} 2> {log} "



## This rule takes the vcf files for each small chunk (scatter) of the chroms and scaffold groups
# and concats them back together based on chrom or scaffold group. So we end up with one vcf file
# per chrom or scaffold group that contains all the samples variant info. 
rule gather_scattered_vcfs:
    input:
        vcf=expand("results/calling/vcf_sections/{sgc}/{scat}.vcf.gz", sgc=sg_or_chrom, scat=unique_scats),
        tbi=expand("results/calling/vcf_sections/{sgc}/{scat}.vcf.gz.tbi", sgc=sg_or_chrom, scat=unique_scats)
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