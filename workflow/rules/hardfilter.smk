# this is the hard-filtering routine recommended by
# the folks that develop the GATK.

# It marks sites as PASS in the FILTER column,
# or it writes out the name of the particular filters
# that were not passed.

rule make_snp_vcf:
    input:
        vcf="results/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz",
        tbi="results/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/snps-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/snps-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/gatk/selectvariants/select-snps-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/make_snp_vcf/selectvariants-snps-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type SNP -O {output.vcf} > {log} 2>&1 "





rule make_indel_vcf:
    input:
        vcf="results/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz",
        tbi="results/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/indels-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/indels-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/gatk/selectvariants/select-indels-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/make_indel_vcf/selectvariants-indels-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type INDEL -O {output.vcf} > {log} 2>&1 "





# filter the snps and indels according to the recommendations at:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
rule hard_filter_snps:
    input:
        vcf="results/hard_filtering/snps-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/snps-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/gatk/variantfiltration/snps-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filter_snps/variantfiltration-snps-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    params:
        filters=config["params"]["gatk"]["VariantFiltration"]["snps"]
    shell:
        "gatk VariantFiltration "
        " -V {input.vcf} "
        " {params.filters} "
        " -O {output.vcf} > {log} 2>&1 "





rule hard_filter_indels:
    input:
        vcf="results/hard_filtering/indels-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/indels-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/gatk/variantfiltration/indels-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filter_indels/variantfiltration-indels-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    params:
        filters=config["params"]["gatk"]["VariantFiltration"]["indels"]
    shell:
        "gatk VariantFiltration "
        " -V {input.vcf} "
        " {params.filters} "
        " -O {output.vcf} > {log} 2>&1 "





rule bung_filtered_vcfs_back_together:
    input:
        snp="results/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz",
        indel="results/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz",
        snp_idx="results/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz.tbi",
        indel_idx="results/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/both-filtered-{sg_or_chrom}.bcf",
    log:
        "results/logs/bung_filtered_vcfs_back_together/bung-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/bung_filtered_vcfs_back_together/bcftools-{sg_or_chrom}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools concat -a {input.snp} {input.indel} | "
        " bcftools view -Ob > {output.vcf}; ) 2> {log} "







rule maf_filter:
    input:
        "results/hard_filtering/both-filtered-{sg_or_chrom}.bcf"
    output:
        "results/hard_filtering/both-filtered-{sg_or_chrom}-maf-{maf}.bcf"
    log:
        "results/logs/maf_filter/{sg_or_chrom}-maf-{maf}.log",
    params:
        maf="{maf}"
    benchmark:
        "results/benchmarks/maf_filter/{sg_or_chrom}-maf-{maf}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools view -Ob -i 'FILTER=\"PASS\" & MAF > {params.maf} ' "
        " {input} > {output} 2>{log} "






### from crct-snake-proj

rule concat_vcfs:
  input:
    vcfs=expand("results/missing-corrected/{c}.vcf.gz", c=CHROMOS)
  output:
    vcf="results/vcf/all.vcf.gz"
  conda:
    "envs/bcftools.yaml"
  log:
    "results/concat_vcfs/all.log"
  shell:
    "bcftools concat -n {input.vcfs} > {output.vcf} 2> {log} "




rule merge_vcfs:
  input:
    indels="results/hard_filtering/indels-filtered/{chromo}.vcf.gz",
    snps="results/hard_filtering/snps-filtered/{chromo}.vcf.gz"
  output:
    vcf="results/hard_filtering/merged/{chromo}.vcf.gz"
  conda:
    "envs/gatk.yaml"
  log:
    "results/merge_vcfs/{chromo}.log"
  benchmark:
    "benchmarks/merge_vcfs/{chromo}.tsv"
  shell:
    "gatk MergeVcfs -I {input.indels} -I {input.snps} -O {output.vcf} 2> {log} "




rule correct_merged_vcfs:
  input:
    vcfs="results/hard_filtering/merged/{chromo}.vcf.gz"
  output:
    vcfs="results/missing-corrected/{chromo}.vcf.gz"
  conda:
    "envs/bcftools.yaml"
  log:
    "results/correct_merged_vcfs/{chromo}.log"
  shell:
    "(bcftools +setGT {input.vcfs} -- -t q -n . -i 'FMT/DP=0 | (FMT/PL[:0]=0 & FMT/PL[:1]=0 & FMT/PL[:2]=0)' | "
    "bcftools +fill-tags - -- -t 'NMISS=N_MISSING' | "
    "bcftools view -Oz - > {output.vcfs}; "
    "bcftools index -t {output.vcfs}) 2> {log} "