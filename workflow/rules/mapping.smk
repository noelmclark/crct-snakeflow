## this rule file includes steps for trimming raw fastqs and mapping them to an indexed reference genome 
# before marking the duplicate reads

rule trim_reads:
    input:
        unpack(get_fastq),
    output:
        #r1=temp("results/trimmed/{sample}---{unit}_R1.fastq.gz"),
        #r2=temp("results/trimmed/{sample}---{unit}_R2.fastq.gz"),
        r1="results/trimmed/{sample}---{unit}_R1.fastq.gz",
        r2="results/trimmed/{sample}---{unit}_R2.fastq.gz",
        html="results/qc/fastp/{sample}---{unit}.html",
        json="results/qc/fastp/{sample}---{unit}.json",
    conda:
        "../envs/fastp.yaml"
    log:
        out="results/logs/trim_reads/{sample}---{unit}.log",
        err="results/logs/trim_reads/{sample}---{unit}.err"
    resources:
        mem_mb=7480,
        time="06:00:00"
    benchmark:
        "results/benchmarks/trim_reads_pe/{sample}---{unit}.bmk"
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
        r1="results/trimmed/{sample}---{unit}_R1.fastq.gz",
        r2="results/trimmed/{sample}---{unit}_R2.fastq.gz",
        genome="resources/genome/OmykA.fasta",
        idx=rules.bwa_index.output,
        #idx=multiext("resources/genome/OmykA.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    output:
        "results/mapped/{sample}---{unit}.sorted.bam"
    conda:
        "envs/bwa2sam.yaml"
    log:
        "results/logs/map_reads/{sample}---{unit}.log"
    benchmark:
        "results/benchmarks/map_reads/{sample}---{unit}.bmk"
    threads: 4
    resources:
        mem_mb=19200,
        time="23:59:59"
    params:
        RG=get_read_group
    shell:
        " (bwa-mem2 mem -t {threads} {params.RG} {input.genome} {input.r1} {input.r2} | "
        " samtools view -u | "
        " samtools sort - > {output}) 2> {log} "


## right now the params in for markduplicates in the config have TAGGING_POLICY All, so duplicates are tagged but not removed
rule mark_duplicates:
    input:
        get_all_bams_of_common_sample
    output:
        bam="results/mkdup/{sample}.bam",
        bai="results/mkdup/{sample}.bai",
        metrics="results/qc/mkdup/{sample}.metrics.txt",
    conda:
        "envs/gatk.yaml"
    log:
        "results/logs/mark_duplicates/{sample}.log",
    benchmark:
        "results/benchmarks/mark_duplicates/{sample}.bmk"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"],
    resources:
        cpus = 1
        mem_mb=112200,
    threads: 30,
    shell:
        " BAMS=$(echo {input} | awk '{{for(i=1;i<=NF;i++) printf(\"-I %s \", $i)}}'); "
        " gatk --java-options '-Xmx3740M' MarkDuplicates  "
        "  {params.extra} "
        "  $BAMS "
        "  -O {output.bam} "
        "  -M {output.metrics} > {log} 2>&1 "