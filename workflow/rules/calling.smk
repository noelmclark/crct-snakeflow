### copied from Eric's calling.smk

# this is primarily a calling thing so we will define it here:
rule make_scaff_group_interval_lists:
    input:
        scaff_groups = config["scaffold_groups"]
    log:
        "results/logs/make_scaff_group_interval_lists/{scaff_group}.log"
    output:
        "results/interval_lists/{scaff_group}.list"
    shell:
        " awk -v sg={wildcards.scaff_group} 'NR>1 && $1 == sg {{print $2}}' {input.scaff_groups} > {output} 2> {log};"


# this is primarily a calling thing so we will define it here:
rule make_chromo_interval_lists:
    log:
        "results/logs/make_chromo_interval_lists/{chromo}.log"
    output:
        "results/interval_lists/{chromo}.list"
    shell:
        " echo {wildcards.chromo} > {output} 2> {log};"


rule make_gvcf_sections:
    input:
        unpack(get_bams_for_calling),
        #bam="results/mkdup/{sample}.bam",
        #bai="results/mkdup/{sample}.bai",
        ref="resources/genome/OmykA.fasta",
        idx="resources/genome/OmykA.dict",
        fai="resources/genome/OmykA.fasta.fai",
        interval_list="results/interval_lists/{sg_or_chrom}.list"
    output:
        gvcf="results/make_gvcf_sections/{sample}/{sg_or_chrom}.g.vcf.gz",
        idx="results/make_gvcf_sections/{sample}/{sg_or_chrom}.g.vcf.gz.tbi",
    conda:
        "../envs/gatk.yaml"
    log:
        stderr="results/logs/make_gvcf_sections/{sample}/{sg_or_chrom}.stderr",
        stdout="results/logs/make_gvcf_sections/{sample}/{sg_or_chrom}.stdout",
    benchmark:
        "results/benchmarks/make_gvcfs_sections/{sample}/{sg_or_chrom}.bmk"
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





rule concat_make_gvcf_sections:
    input: 
        expand("results/make_gvcf_sections/{{sample}}/{sgc}.g.vcf.gz", sgc = unique_chromosomes + unique_scaff_groups)
    output:
        gvcf="results/gvcf/{sample}.g.vcf.gz",
        idx="results/gvcf/{sample}.g.vcf.gz.tbi"
    log:
        "results/logs/concat_make_gvcf_sections/{sample}.txt"
    benchmark:
        "results/benchmarks/concat_make_gvcf_sections/{sample}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools concat {params.opts} -O z {input} > {output.gvcf} 2>{log}; "
        " bcftools index -t {output.gvcf} "






rule import_genomics_db_by_chromo:
    input:
        gvcfs=expand("results/gvcf/{{chromo}}/{s}.g.vcf.gz", s=SAMPLES)
    output:
        gdb=directory("results/genomics_db/{chromo}")
    conda:
        "../envs/gatk.yaml"
    log:
        "results/logs/import_genomics_db_by_chromo/{chromo}.log"
    benchmark:

    params:
        java_opts="-Xmx4g"
    shell:
        " VS=$(for i in {input.gvcfs}; do echo -V $i; done); "  # make a string like -V file1 -V file2
        " gatk --java-options \"-Xmx4g\" GenomicsDBImport "
        "  $VS  "
        "  --genomicsdb-workspace-path {output.gdb} "
        "  -L  {wildcards.chromo} 2> {log} "

# apparently setting Xmx4G affects the resource allocation here, but
# I want to be sure to get more RAM than that.  On slurm, it will get that
# if we request 2 CPUs on sedna.  Note that this thing runs no faster with
# four cpus and reader threads than with 2.
## Note: the next two rules are unsatisfying because params.fileflags and input.gvcfs are
## defined separately, though they are effecively the same thing.
## This is pretty messed up.  The log file is for the last run, but I also
## tee it to one named with the import number.  
rule genomics_db_import_chromosomes:
    input:
        gvcfs=lambda wc: expand("results/make_gvcf_sections/{sample}/{{chromo}}.g.vcf.gz", sample=sample_list),
        gvcf_idxs=lambda wc: expand("results/make_gvcf_sections/{sample}/{{chromo}}.g.vcf.gz.tbi", sample=sample_list),
    output:
        db=directory("results/genomics_db/{chromo}")
    log:
        "results/logs/gatk/genomicsdbimport/{chromo}.log"
    benchmark:
        "results/benchmarks/genomics_db_import_chromosomes/{chromo}.bmk"
    params:
        my_opts=chromo_import_gdb_opts,
        java_opts="-Xmx4g",  # optional
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk --java-options {params.java_opts} GenomicsDBImport "
        " $(echo {input.gvcfs} | awk '{{for(i=1;i<=NF;i++) printf(\" -V %s \", $i)}}') "
        " {params.my_opts} {output.db} > {log} 2>&1  "
        


# It seems like I have to set threads to 2 in order to get 2 CPUs for this
# job.  resources: cpus = 2, does not work.  Same for above. I set reader
# threads to 2.  It doesn't seem that there is any advantage to doing more
# than that.
rule genomics_db_import_scaffold_groups:
    input:
        gvcfs=lambda wc: expand("results/make_gvcf_sections/{sample}/{{scaff_group}}.g.vcf.gz", sample=sample_list),
        gvcf_idxs=lambda wc: expand("results/make_gvcf_sections/{sample}/{{scaff_group}}.g.vcf.gz.tbi", sample=sample_list),
        scaff_groups = config["scaffold_groups"],
    output:
        interval_list="results/gdb_intervals/{scaff_group}.list",
        db=directory("results/genomics_db/{scaff_group}")
    log:
        "results/logs/gatk/genomicsdbimport/{scaff_group}.log"
    benchmark:
        "results/benchmarks/genomics_db_import_scaffold_groups/{scaff_group}.bmk"
    params:
        my_opts=scaff_group_import_gdb_opts,
        java_opts="-Xmx4g",  # optional
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " awk -v sg={wildcards.scaff_group} 'NR>1 && $1 == sg {{print $2}}' {input.scaff_groups} > {output.interval_list}; "
        " gatk --java-options {params.java_opts} GenomicsDBImport "
        " $(echo {input.gvcfs} | awk '{{for(i=1;i<=NF;i++) printf(\" -V %s \", $i)}}') "
        " {params.my_opts} {output.db} >{log} 2>&1; "





rule vcf_from_gdb_by_chromo:
    input:
        gdb="results/genomics_db/{chromo}",
        ref="resources/genome/OmykA.fasta",
        fai="resources/genome/OmykA.fasta.fai",
        idx="resources/genome/OmykA.dict",
    output:
        vcf="results/chromo_vcfs/{chromo}.vcf.gz",
        idx="results/chromo_vcfs/{chromo}.vcf.gz.tbi",
    conda:
        "../envs/gatk.yaml"
    log:
        "results/logs/vcf_from_gdb_by_chromo/{chromo}.txt"
    shell:
        " gatk --java-options \"-Xmx4g\" GenotypeGVCFs "
        "  -R {input.ref}  "
        "  -V gendb://{input.gdb} "
        "  -O {output.vcf} 2> {log} "





# {type_of_subset} will be either "chromosomes" or "scaffold_groups"
# {sg_or_chrom} will be either like "CM031199.1" (if type_of_subset is chromosomes), 
# or {scaff_group001} if type_of_subset is scaffold_groups.  
rule genomics_db2vcf_scattered:
    input:
        genome="resources/genome.fasta",
        scatters="results/scatter_interval_lists/{sg_or_chrom}/{scatter}.list",
        db="results/genomics_db/{sg_or_chrom}",
    output:
        vcf="results/vcf_sections/{sg_or_chrom}/{scatter}.vcf.gz",
        tbi="results/vcf_sections/{sg_or_chrom}/{scatter}.vcf.gz.tbi"
    log:
        "results/logs/gatk/genotypegvcfs/{sg_or_chrom}/{scatter}.log",
    benchmark:
        "results/benchmarks/genomics_db2vcf/{sg_or_chrom}/{scatter}.bmk",
    params:
        gendb="results/genomics_db/{sg_or_chrom}",
        java_opts="-Xmx8g",  # I might need to consider a temp directory, too in which case, put it in the config.yaml
        pextra=" --genomicsdb-shared-posixfs-optimizations --only-output-calls-starting-in-intervals "
    resources:
        mem_mb = 11750,
        cpus = 2,
        time = "1-00:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk --java-options {params.java_opts} GenotypeGVCFs "
        " {params.pextra} "
        " -L {input.scatters} "
        " -R {input.genome} "
        " -V gendb://{params.gendb} "
        " -O {output.vcf} > {log} 2> {log} "


rule gather_scattered_vcfs:
    input:
        vcf=lambda wc: get_scattered_vcfs(wc, ""),
        tbi=lambda wc: get_scattered_vcfs(wc, ".tbi"),
    output:
        vcf="results/vcf_sections/{sg_or_chrom}.vcf.gz",
        tbi="results/vcf_sections/{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/gather_scattered_vcfs/{sg_or_chrom}.txt"
    benchmark:
        "results/benchmarks/gather_scattered_vcfs/{sg_or_chrom}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Oz {input.vcf} > {output.vcf}; "
        " bcftools index -t {output.vcf})  2>{log}; "




# this is a little rule we throw in here so that we can mark
# an individual as missing data (./. or .|.) when it has a read
# depth of 0, because GATK now marks those as 0/0,
# see https://gatk.broadinstitute.org/hc/en-us/community/posts/4476803114779-GenotypeGVCFs-Output-no-call-as-reference-genotypes?page=1#community_comment_6006727219867
# this also adds an INFO field NMISS, which gives the number of samples missing a call.
# 8/26/22: This has been updated to also mark genotypes as missing if they have a PL of 0,0,0.
rule mark_dp0_as_missing:
    input:
        vcf="results/vcf_sections/{sg_or_chrom}.vcf.gz"
    output:
        vcf="results/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz",
        tbi="results/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/mark_dp0_as_missing/{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/mark_dp0_as_missing/{sg_or_chrom}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools +setGT {input.vcf} -- -t q -n . -i 'FMT/DP=0 | (FMT/PL[:0]=0 & FMT/PL[:1]=0 & FMT/PL[:2]=0)' | "
        " bcftools +fill-tags - -- -t 'NMISS=N_MISSING' | "
        " bcftools view -Oz - > {output.vcf}; "
        " bcftools index -t {output.vcf}) 2> {log} "




rule bcf_concat:
    input:
        expand("results/hard_filtering/both-filtered-{sgc}.bcf", sgc = unique_chromosomes + unique_scaff_groups)
    output:
        bcf="results/bcf/all.bcf",
        tbi="results/bcf/all.bcf.csi"
    log:
        "results/logs/bcf_concat/bcf_concat_log.txt"
    benchmark:
        "results/benchmarks/bcf_concat/bcf_concat.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2> {log}; "





rule bcf_concat_mafs:
    input:
        expand("results/hard_filtering/both-filtered-{sgc}-maf-{{maf}}.bcf", sgc = unique_chromosomes + unique_scaff_groups)
    output:
        bcf="results/bcf/pass-maf-{maf}.bcf",
        tbi="results/bcf/pass-maf-{maf}.bcf.csi"
    log:
        "results/logs/bcf_concat_mafs/maf-{maf}.txt"
    benchmark:
        "results/benchmarks/bcf_concat_mafs/maf-{maf}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2>{log}; "