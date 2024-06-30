# this is the hard-filtering routine recommended by
# the folks that develop the GATK.

# It marks sites as PASS in the FILTER column,
# or it writes out the name of the particular filters
# that were not passed.

# selects only snps by chrom or scaffold group from the missing-corrected vcf sections
rule make_snp_vcf:
    input:
        vcf="results/calling/corrected_missing_vcf_sect/{sg_or_chrom}.vcf.gz",
        tbi="results/calling/corrected_missing_vcf_sect/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/snps-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/snps-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/hard_filtering/selectvariants/select-snps-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/selectvariants/select-snps-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type SNP -O {output.vcf} > {log} 2>&1 "





#select only indels by chrom or scaffold group from the missing-corrected vcf sections
rule make_indel_vcf:
    input:
        vcf="results/calling/corrected_missing_vcf_sect/{sg_or_chrom}.vcf.gz",
        tbi="results/calling/corrected_missing_vcf_sect/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/indels-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/indels-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/hard_filtering/selectvariants/select-indels-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/selectvariants/select-indels-{sg_or_chrom}.bmk"
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
        "results/logs/hard_filtering/variantfiltration/hard-filter-snps-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/variantfiltration/hard-filter-snps-{sg_or_chrom}.bmk"
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
        "results/logs/hard_filtering/variantfiltration/hard-filter-indels-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/variantfiltration/hard-filter-indels-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    params:
        filters=config["params"]["gatk"]["VariantFiltration"]["indels"],
    shell:
        "gatk VariantFiltration "
        " -V {input.vcf} "
        " {params.filters} "
        " -O {output.vcf} > {log} 2>&1 "



# merge the snp and indel vcfs that were hard filtered back together by chrom or scaffold_groups
rule merge_filtered_vcfs:
    input:
        snp="results/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz",
        indel="results/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz",
        snp_idx="results/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz.tbi",
        indel_idx="results/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/hardfiltered-merged-{sg_or_chrom}.vcf.gz",
    log:
        "results/logs/hard_filtering/merge_filtered_vcfs/hardfiltered-merged-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/merge_filtered_vcfs/hardfiltered-merged-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    shell:
        " gatk MergeVcfs -I {input.indel} -I {input.snp} -O {output.vcf} 2> {log} "
  


# this gives us a single merged bcf file without any maf filtering
rule bcf_concat_hardfilter_merged_sect:
    input:
        expand("results/hard_filtering/hardfiltered-merged-{sgc}.vcf.gz", sgc=sg_or_chrom),
    output:
        bcf="results/bcf/all.bcf",
        tbi="results/bcf/all.bcf.csi"
    log:
        "results/logs/bcf_concat_hardfilter_merged_sect/bcf_concat_log.txt"
    benchmark:
        "results/benchmarks/bcf_concat_hardfilter_merged_sect/bcf_concat.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2> {log}; "



## this will run on each scaffold group or chromosome 
# for each of the maf options specified in the maf_cutoff part of the config.yaml file
# I have it set up to run on the missing-corrected-merged-vcf files from the correct_merged_vcfs rule 
rule maf_filter:
    input:
        "results/hard_filtering/hardfiltered-merged-{sg_or_chrom}.vcf.gz"
    output:
        "results/hard_filtering/hardfiltered-merged-{sg_or_chrom}-maf-{mafs}.bcf"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/maf_filter/merged-{sg_or_chrom}-maf-{mafs}.log",
    benchmark:
        "results/benchmarks/hard_filtering/maf_filter/merged-{sg_or_chrom}-maf-{mafs}.bmk"
    params:
        maf="{mafs}"
    shell:
        " bcftools view -Ob -i 'FILTER=\"PASS\" & MAF > {params.maf} ' "
        " {input} > {output} 2>{log} "



## this gives us a single merged bcf file for each of the maf cuttoffs we specify that contains only variants 
# which pass the specified maf cuttoff
rule concat_mafs_bcf:
    input:
        expand("results/hard_filtering/hardfiltered-merged-{sgc}-maf-{{maf}}.bcf", sgc=sg_or_chrom)
    output:
        bcf="results/bcf/all-that-pass-maf-{maf}.bcf",
        tbi="results/bcf/all-that-pass-maf-{maf}.bcf.csi"
    log:
        "results/logs/hard_filtering/concat_mafs_bcf/maf-{maf}.txt"
    benchmark:
        "results/benchmarks/hard_filtering/concat_mafs_bcf/maf-{maf}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2>{log}; "