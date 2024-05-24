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






# eca modified this.  The idea is to give 4 threads to bwa.
# and it will get 4 cores and also take all the memory you'd
# expect for those cores.  Sedna's machines are almost all
# 20 core units, so this should fill them up OK.
rule map_reads:
    input:
        reads = [
            "results/trimmed/{sample}---{unit}.1.fastq.gz",
            "results/trimmed/{sample}---{unit}.2.fastq.gz"
        ],
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/{sample}---{unit}.sorted.bam"),
    log:
        "results/logs/map_reads/{sample}---{unit}.log",
    benchmark:
        "results/benchmarks/map_reads/{sample}---{unit}.bmk"
    params:
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
        sort_extra=""
    threads: 4
    resources:
        mem_mb=19200,
        time="23:59:59"
    wrapper:
        "v1.23.3/bio/bwa/mem"



rule mark_duplicates:
    input:
        get_all_bams_of_common_sample
    output:
        bam="results/mkdup/{sample}.bam",
        bai="results/mkdup/{sample}.bai",
        metrics="results/qc/mkdup/{sample}.metrics.txt",
    log:
        "results/logs/picard/mkdup/{sample}.log",
    benchmark:
        "results/benchmarks/mark_duplicates/{sample}.bmk"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"],
    resources:
        cpus = 1
    wrapper:
        "v1.1.0/bio/picard/markduplicates"








rule map_reads:
  input:
    r1="results/trimmed/{sample}---{unit}_R1.fastq.gz",
    r2="results/trimmed/{sample}---{unit}_R2.fastq.gz",
    genome="data/genome/OmykA.fasta",
    idx=multiext("data/genome/OmykA.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
  output:
    "results/mapped/{sample}---{unit}.sorted.bam"
  conda:
    "envs/bwa2sam.yaml"
  log:
    "results/logs/map_reads/{sample}---{unit}.log"
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



rule mark_duplicates:
  input:
    bam=get_all_bams_of_common_sample,
  output:
    bam="results/mkdup/{sample}.bam",
    bai="results/mkdup/{sample}.bai",
    metrics="results/qc/mkdup_metrics/{sample}.metrics"
  conda:
    "envs/gatk.yaml"
  log:
    "results/logs/mark_duplicates/{sample}.log"
  resources:
    mem_mb=112200,
  threads: 30,
  #params:
  #  extra=config["params"]["picard"]["MarkDuplicates"],
  #wrapper:
  #  "v3.9.0/bio/picard/markduplicates"
  shell:
    " BAMS=$(echo {input.bam} | awk '{{for(i=1;i<=NF;i++) printf(\"-I %s \", $i)}}'); "
    " gatk --java-options '-Xmx3740M' MarkDuplicates  "
    "  --CREATE_INDEX --TMP_DIR results/snake-tmp "
    "  $BAMS "
    "  -O {output.bam} "
    "  -M {output.metrics} > {log} 2>&1 "