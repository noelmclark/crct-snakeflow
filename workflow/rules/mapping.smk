## this rule file includes steps for trimming raw fastqs and mapping them to an indexed reference genome 
# before marking the duplicate reads

rule trim_reads:
    input:
        unpack(get_fastqs),
    output:
        r1=temp("results/mapping/trimmed/{sample}---{unit}_R1.fastq.gz"),
        r2=temp("results/mapping/trimmed/{sample}---{unit}_R2.fastq.gz"),
        html="results/qc/fastp/{sample}---{unit}.html",
        json="results/qc/fastp/{sample}---{unit}.json",
    conda:
        "../envs/fastp.yaml"
    log:
        out="results/logs/mapping/trim_reads/{sample}---{unit}.log",
        err="results/logs/mapping/trim_reads/{sample}---{unit}.err"
    resources:
        mem_mb=7480,
        time="06:00:00"
    benchmark:
        "results/benchmarks/mapping/trim_reads/{sample}---{unit}.bmk"
    params:
        as1=config["params"]["fastp"]["adapter_sequence1"],
        as2=config["params"]["fastp"]["adapter_sequence2"],
        parm=config["params"]["fastp"]["other_options"]
    shell:
        " fastp -i {input.r1} -I {input.r2}       "
        "       -o {output.r1} -O {output.r2}     "
        "       -h {output.html} -j {output.json} "
        "  --adapter_sequence={params.as1}        "
        "  --adapter_sequence_r2={params.as2}     "
        "  {params.parm} > {log.out} 2> {log.err} "



rule map_reads:
    input:
        r1="results/mapping/trimmed/{sample}---{unit}_R1.fastq.gz",
        r2="results/mapping/trimmed/{sample}---{unit}_R2.fastq.gz",
        genome="resources/genome/OmykA.fasta",
        idx=multiext("resources/genome/OmykA.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        temp("results/mapping/mapped/{sample}---{unit}.sorted.bam"),
    conda:
        "../envs/bwa2sam.yaml"
    log:
        "results/logs/mapping/map_reads/{sample}---{unit}.log",
    benchmark:
        "results/benchmarks/mapping/map_reads/{sample}---{unit}.bmk",
    threads: 8
    resources:
        mem_mb=30000,
        time="23:59:59",
    params:
        RG=get_read_group
    shell:
        " (bwa-mem2 mem -t {threads} {params.RG} {input.genome} {input.r1} {input.r2} | "
        " samtools view -u | "
        " samtools sort - > {output}) 2> {log} "


## right now the params in for markduplicates in the config have TAGGING_POLICY All, so duplicates are tagged but not removed
# if you want to remove duplicates, use the rule below this one
# also, if you assigned your read groups properly, input should be get_all_bams_of_common_sample instead...
rule mark_duplicates:
    input:
        get_all_bams_of_common_sample #or get_RG_fixed_bams_of_common_sample
    output:
        bam="results/mapping/mkdup/mkdup-{sample}.bam",
        bai="results/mapping/mkdup/mkdup-{sample}.bai",
        metrics="results/qc/mkdup/mkdup-{sample}.metrics.txt",
    conda:
        "../envs/gatk.yaml"
    log:
        "results/logs/mapping/mark_duplicates/{sample}.log",
    benchmark:
        "results/benchmarks/mapping/mark_duplicates/{sample}.bmk"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"],
    resources:
        cpus = 1,
        mem_mb=112200,
    threads: 30,
    shell:
        " BAMS=$(echo {input} | awk '{{for(i=1;i<=NF;i++) printf(\"-I %s \", $i)}}'); "
        " gatk --java-options '-Xmx3740M' MarkDuplicates  "
        "  {params.extra} "
        "  $BAMS "
        "  -O {output.bam} "
        "  -M {output.metrics} > {log} 2>&1 "

## this rule merges bams of common sample and removes the duplicates (--REMOVE_DUPLICATES true) 
# also, if you assigned your read groups properly, input should be get_all_bams_of_common_sample instead...
rule GATK_remove_duplicates:
    input:
        get_all_bams_of_common_sample #or get_RG_fixed_bams_of_common_sample
    output:
        bam="results/mapping/gatk-rmdup/{sample}.bam",
        bai="results/mapping/gatk-rmdup/{sample}.bai",
        metrics="results/qc/gatk-rmdup/{sample}.metrics.txt",
    conda:
        "../envs/gatk.yaml"
    log:
        "results/logs/mapping/GATK_remove_duplicates/{sample}.log",
    benchmark:
        "results/benchmarks/mapping/GATK_remove_duplicates/{sample}.bmk"
    params:
        "--REMOVE_DUPLICATES true --CREATE_INDEX --TMP_DIR results/snake-tmp "
    resources:
        cpus = 1,
        mem_mb=112200,
    threads: 30,
    shell:
        " BAMS=$(echo {input} | awk '{{for(i=1;i<=NF;i++) printf(\"-I %s \", $i)}}'); "
        " gatk --java-options '-Xmx3740M' MarkDuplicates  "
        "  {params} "
        "  $BAMS "
        "  -O {output.bam} "
        "  -M {output.metrics} > {log} 2>&1 "

## rule to remove the duplicates from the marked_duplicates bams to be used for downstream analyses
# an alternative to removing them with GATK MarkDuplicates
# -F 1024 excludes reads with flag 1024 (optical or PCR duplicate)
# -F 3852 was used in the Leopard paper 
# https://broadinstitute.github.io/picard/explain-flags.html 
#rule samtools_remove_duplicates:
#    input:
#        bam="results/mapping/mkdup/mkdup-{sample}.bam",
#        bai="results/mapping/mkdup/mkdup-{sample}.bai",
#    output:
#        bam="results/mapping/samtools-rmdup/{sample}.bam",
#        bai="results/mapping/samtools-rmdup/{sample}.bai",
#    log:
#        rm="results/logs/mapping/samtools-rmdup/{sample}.log",
#        index="results/logs/mapping/samtools-rmdup/index-{sample}.log",
#    conda:
#        "../envs/bamutil_samtools.yaml"
#    benchmark:
#        "results/benchmarks/mapping/samtools-rmdup/{sample}.bmk",
#    shell:
#        " samtools view -F 1024 {input.bam} -o {output.bam} >2 {log.rm} && "
#        " samtools index {output.bai} 2> {log.index} "
