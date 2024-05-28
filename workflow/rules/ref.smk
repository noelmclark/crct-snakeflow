# Rule to download the genome using wget from NCBI if you don't have the fasta locally  
#rule get_genome:
#    output:
#        "resources/genome/OmykA.fasta",
#    log:
#        "results/logs/get_OmykA.log",
#    benchmark:
#        "results/benchmarks/get_genome/get_OmykA.bmk",
#    params:
#        url=config["ref"]["genome_url"],
#    conda:
#        "../envs/wget.yaml"
#    shell:
#        " (tmp_dir=$(mktemp -d) && "
#        " URL={params.url} && "
#        " if [[ $URL =~ \.gz$ ]]; then EXT='.gz'; else EXT=''; fi && "
#        " wget -O $tmp_dir/file$EXT $URL && "
#        " if [[ $URL =~ \.gz$ ]]; then gunzip $tmp_dir/file$EXT; fi && "
#        " mv $tmp_dir/file {output}) > {log} 2>&1 "



rule genome_faidx:
    input:
        "resources/genome/OmykA.fasta",
    output:
        "resources/genome/OmykA.fasta.fai",
    log:
        "results/logs/genome_faidx.log",
    benchmark:
        "results/benchmarks/genome_faidx/genome_faidx.bmk",
    conda:
        "../envs/bwa2sam.yaml"
    shell:
        "samtools faidx {input} 2> {log} "



rule genome_dict:
    input:
        "resources/genome/OmykA.fasta",
    output:
        "resources/genome/OmykA.dict",
    log:
        "results/logs/genome_dict.log",
    benchmark:
        "results/benchmarks/genome_dict/genome_dict.bmk"
    conda:
        "../envs/bwa2sam.yaml"
    shell:
        "samtools dict {input} > {output} 2> {log} "



rule bwa_index:
    input:
        "resources/genome/OmykA.fasta",
    output:
        multiext("resources/genome/OmykA.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        out="results/logs/bwa_index.log",
        err="results/logs/bwa_index.err",
    benchmark:
        "results/benchmarks/bwa_index/bwa_index.bmk",
    resources:
        mem_mb=75800,
        time="12:00:00"
    conda:
        "../envs/bwa2sam.yaml"
    shell:
        " bwa-mem2 index {input} > {log.out} 2> {log.err} "
