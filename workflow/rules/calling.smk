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
        bam="results/mapping/mkdup/mkdup-{sample}.bam",
        bai="results/mapping/mkdup/mkdup-{sample}.bai",
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
        expand("results/calling/gvcf_sections/{s}/{sgc}.g.vcf.gz", s=sample, sgc=sg_or_chrom)
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
