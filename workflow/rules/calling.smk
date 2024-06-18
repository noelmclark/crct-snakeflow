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
        bam="results/angsd_bams/overlap_clipped/{sample}.bam",
        bai="results/angsd_bams/overlap_clipped/{sample}.bai",
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
        java_opts="-Xmx4g"
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    shell:
        " VS=$(for i in {input.gvcfs}; do echo -V $i; done); "  # make a string like -V file1 -V file2
        " gatk --java-options {params.java_opts} GenomicsDBImport "
        "  $VS  "
        #" $(echo {input.gvcfs} | awk '{{for(i=1;i<=NF;i++) printf(\" -V %s \", $i)}}') "
        "  --genomicsdb-workspace-path {output.gdb} "
        "  -L {wildcards.chromo} 2> {log} "


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
        java_opts="-Xmx4g"
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    shell:
        " VS=$(for i in {input.gvcfs}; do echo -V $i; done); "  # make a string like -V file1 -V file2
        " gatk --java-options {params.java_opts} GenomicsDBImport "
        "  $VS  "
        "  --genomicsdb-workspace-path {output.gdb} "
        "  -L {wildcards.scaff_group} 2> {log} "


