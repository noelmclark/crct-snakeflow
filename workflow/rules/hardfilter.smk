#### This rule file contains many rule chunks related to filtering at different levels ####

### these first 6 rules are the hard-filtering routine recommended by the folks that develop the GATK. ###
# It marks sites as PASS in the FILTER column,
# or it writes out the name of the particular filters
# that were not passed.

# selects only snps by chrom or scaffold group from the missing-corrected vcf sections
rule make_snp_vcf:
    input:
        vcf="results/calling/corrected_missing_vcf_sect/{sg_or_chrom}.vcf.gz",
        tbi="results/calling/corrected_missing_vcf_sect/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/select-variants/snps-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/select-variants/snps-{sg_or_chrom}.vcf.gz.tbi"
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
        vcf="results/hard_filtering/select-variants/indels-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/select-variants/indels-{sg_or_chrom}.vcf.gz.tbi"
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
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
rule hard_filter_snps:
    input:
        vcf="results/hard_filtering/select-variants/snps-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/select-variants/snps-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/variant-filtered/snps-filtered-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/variant-filtered/snps-filtered-{sg_or_chrom}.vcf.gz.tbi"
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
        vcf="results/hard_filtering/select-variants/indels-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/select-variants/indels-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/variant-filtered/indels-filtered-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/variant-filtered/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
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
        snp="results/hard_filtering/variant-filtered/snps-filtered-{sg_or_chrom}.vcf.gz",
        indel="results/hard_filtering/variant-filtered/indels-filtered-{sg_or_chrom}.vcf.gz",
        snp_idx="results/hard_filtering/variant-filtered/snps-filtered-{sg_or_chrom}.vcf.gz.tbi",
        indel_idx="results/hard_filtering/variant-filtered/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/merged-filtered/hardfiltered-merged-{sg_or_chrom}.vcf.gz",
    log:
        "results/logs/hard_filtering/merge_filtered_vcfs/hardfiltered-merged-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/merge_filtered_vcfs/hardfiltered-merged-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    shell:
        " gatk MergeVcfs -I {input.indel} -I {input.snp} -O {output.vcf} 2> {log} "


## convert each merged filtered vcf into a bcf
# I think there's a way to do this all in one with concat but I was running into errors
rule merged_filtered_vcfs_to_bcfs:
    input:
        vcf="results/hard_filtering/merged-filtered/hardfiltered-merged-{sg_or_chrom}.vcf.gz",
    output:
        bcf="results/hard_filtering/merged-filtered-bcf/hardfiltered-merged-{sg_or_chrom}.bcf",
    log:
        "results/logs/hard_filtering/merged_filtered_vcfs_to_bcfs/hardfiltered-merged-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/merged_filtered_vcfs_to_bcfs/hardfiltered-merged-{sg_or_chrom}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools view -Ob {input} > {output} 2>{log} "

##########################################
#### end of GATK hard filtering rules ####
##########################################


## this gives us a single merged bcf file without any filtering & with variants that failed the hard filtering flagged BUT NOT REMVOED
rule bcf_concat_hardfilter_merged_sect:
    input:
        expand("results/hard_filtering/merged-filtered-bcf/hardfiltered-merged-{sgc}.bcf", sgc=sg_or_chrom),
    output:
        bcf="results/bcf/all.bcf",
        tbi="results/bcf/all.bcf.csi"
    log:
        "results/logs/hard_filtering/bcf_concat_hardfilter_merged_sect/bcf_concat_log.txt"
    benchmark:
        "results/benchmarks/hard_filtering/bcf_concat_hardfilter_merged_sect/bcf_concat.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2> {log}; "

###############################


### these next 2 rules give us a single merged bcf file with variants that failed the GATK hard filtering REMVOED ###
# first filter out the variants that failed GATK filters
rule filter_all:
    input:
        "results/hard_filtering/merged-filtered/hardfiltered-merged-{sg_or_chrom}.vcf.gz"
    output:
        "results/hard_filtering/filter-pass/{sg_or_chrom}.bcf"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/filter_pass/{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/filter_pass/{sg_or_chrom}.bmk"
    shell:
        " bcftools view -Ob -i 'FILTER=\"PASS\"' "
        " {input} > {output} 2> {log} "

# then concat all these filtered sections together 
rule concat_filter_bcf:
    input:
        expand("results/hard_filtering/filter-pass/{sgc}.bcf", sgc=sg_or_chrom)
    output:
        bcf="results/bcf/all-filtered.bcf",
        tbi="results/bcf/all-filtered.bcf.csi"
    log:
        "results/logs/hard_filtering/concat_bisnp_bcf/all-filtered.log"
    benchmark:
        "results/benchmarks/hard_filtering/concat_bisnp_bcf/all-filtered.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2>{log}; "

###############################


### these couple rules take the vcf sections and remove all variants that failed the gatk filters ###
## before removing any variants that map to regions outside of autosomes
## they DO NOT perform any MAF filtering -- I will do that in PLINK (11/21/24)

# first filter out the variants that failed GATK filters and include just binsps for each vcf chunk
rule filter_bisnp:
    input:
        "results/hard_filtering/merged-filtered/hardfiltered-merged-{sg_or_chrom}.vcf.gz"
    output:
        "results/hard_filtering/filter-bisnp/bisnps-{sg_or_chrom}.bcf"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/filter_bisnp/bisnps-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/filter_bisnp/bisnps-{sg_or_chrom}.bmk"
    shell:
        " bcftools view -Ob -v snps -m 2 -M 2 -i 'FILTER=\"PASS\"' "
        " {input} > {output} 2> {log} "

# then merge bisnp filtered vcf chunks into bcf file 
rule concat_bisnp_bcf:
    input:
        expand("results/hard_filtering/filter-bisnp/bisnps-{sgc}.bcf", sgc=sg_or_chrom)
    output:
        bcf="results/bcf/bisnps.bcf",
        tbi="results/bcf/bisnps.bcf.csi"
    log:
        "results/logs/hard_filtering/concat_bisnp_bcf/bisnps.log"
    benchmark:
        "results/benchmarks/hard_filtering/concat_bisnp_bcf/bisnps.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2>{log}; "

# last, remove any variants not located on assembled autosomal chroms
rule get_aut_bisnp_bcf:
    input:
        bcf="results/bcf/bisnps.bcf",
        tbi="results/bcf/bisnps.bcf.csi",
        regions="results/psmc/get-aut-regions/autosomal_regions.bed"
    output:
        bcf="results/bcf/aut-bisnps.bcf",
        tbi="results/bcf/aut-bisnps.bcf.csi"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/aut_bisnp_bcf/aut-bisnps.log"
    benchmark:
        "results/benchmarks/hard_filtering/aut_bisnp_bcf/aut-bisnps.bmk"
    shell:
        " (bcftools view -Ob -R {input.regions} {input.bcf} > {output.bcf}; "
        " bcftools index {output.bcf}) 2> {log} "

# this next rule is related, and removes all SNPs within +/- 5bp of indels 
rule filter_snps_near_indels:
    input:
        bcf="results/bcf/aut-bisnps.bcf",
        tbi="results/bcf/aut-bisnps.bcf.csi"
    output:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        tbi="results/bcf/aut-bisnps-no5indel.bcf.csi"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/aut_bisnp_bcf/aut-bisnps-no5indel.log"
    benchmark:
        "results/benchmarks/hard_filtering/aut_bisnp_bcf/aut-bisnps-no5indel.bmk"
    shell:
        " (bcftools filter -Ob --SnpGap 5 {input.bcf} > {output.bcf}; "
        " bcftools index {output.bcf}) 2> {log} "

###############################


### these rules filter for variants that pass the gatk filters and for maf>= 0.05 using bcftools ###

## this will run on each scaffold group or chromosome 
# for each of the maf options specified in the maf_cutoff part of the config.yaml file
# I have it set up to run on the missing-corrected-merged-vcf files from the correct_merged_vcfs rule 
rule filter_bisnp_maf:
    input:
        "results/hard_filtering/merged-filtered/hardfiltered-merged-{sg_or_chrom}.vcf.gz"
    output:
        "results/hard_filtering/filter-bisnp-maf-{maf}/biallelic-snps-{sg_or_chrom}-maf-{mafs}.bcf"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/filter_bisnp_maf_{maf}/biallelic-snps-{sg_or_chrom}-maf-{mafs}.log",
    benchmark:
        "results/benchmarks/hard_filtering/filter_bisnp_maf_{maf}/biallelic-snps-{sg_or_chrom}-maf-{mafs}.bmk"
    params:
        maf="{mafs}"
    shell:
        " bcftools view -Ob -v snps -m 2 -M 2 -i 'FILTER=\"PASS\" && MAF >= 0.05' "
        " {input} > {output} 2>{log} "


## this gives us a single merged bcf file for each of the maf cuttoffs we specify that contains only variants 
# which pass the specified maf cuttoff & all other filters
rule concat_bisnp_maf_bcf:
    input:
        expand("results/hard_filtering/filter-bisnp-maf-{{maf}}/biallelic-snps-{sgc}-maf-{{maf}}.bcf", sgc=sg_or_chrom)
    output:
        bcf="results/bcf/biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/biallelic-snps-maf-{maf}.bcf.csi"
    log:
        "results/logs/hard_filtering/concat_bisnp_mafs_bcf/biallelic-snps-maf-{maf}.log"
    benchmark:
        "results/benchmarks/hard_filtering/concat_bisnp_mafs_bcf/biallelic-snps-maf-{maf}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2>{log}; "


## this rule filters the biallelic snp and maf file to include only autosomal regions and not scaffold reads
rule get_aut_bisnp_mafs_bcf:
    input:
        bcf="results/bcf/biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/biallelic-snps-maf-{maf}.bcf.csi",
        regions="results/psmc/get-aut-regions/autosomal_regions.bed"
    output:
        bcf="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf.csi"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/aut_bisnp_mafs_bcf/autosomal-biallelic-snps-maf-{maf}.log"
    benchmark:
        "results/benchmarks/hard_filtering/aut_bisnp_mafs_bcf/autosomal-biallelic-snps-maf-{maf}.bmk"
    shell:
        " (bcftools view -Ob -R {input.regions} {input.bcf} > {output.bcf}; "
        " bcftools index {output.bcf}) 2> {log} "

###############################
### end of main bcf section ###
###############################



### bcf subset section ###
### These rules filter the bcf down to a desired subset of samples given in a txt file ###

# this part of the workflow was heavily influenced by Eric's mega-post-bcf-exploratory-snakeflows github (https://github.com/eriqande/mega-post-bcf-exploratory-snakeflows)
# it will run on each scaffold group or chromosome for the missing-corrected-merged-vcf files from the correct_merged_vcfs rule 
# including only samples that are listed in the subsamp file 
# the +fill-tags refills all of the tags - this makes it so that our new MAFs in particular are only dependent on the variants present in any of the new subset of samples
# before we filter our variants by a MAF >= 0.05
rule subset_filter_bcf:
    input:
        bcf="results/hard_filtering/merged-filtered/hardfiltered-merged-{sg_or_chrom}.vcf.gz",
        subsamp=get_bcf_subsamp_path
    output:
        temp("results/hard_filtering/subsets/{bcfsubsamp}/bisnps-{sg_or_chrom}-maf-0.05-{bcfsubsamp}.bcf")
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/subsets/{bcfsubsamp}/filter-{sg_or_chrom}-{bcfsubsamp}.log",
    benchmark:
        "results/benchmarks/hard_filtering/subsets/{bcfsubsamp}/filter-{sg_or_chrom}-{bcfsubsamp}.bmk"
    shell:
        " bcftools view -Ou -S {input.subsamp} {input.bcf} | "
        " bcftools +fill-tags -Ou -- -t all | "
        " bcftools view -Ob -v snps -m 2 -M 2 -i 'FILTER=\"PASS\" && MAF >= 0.05' "
        " > {output} 2>{log} "

## this merges the subsamped & filtered bcf file 
rule concat_subset_filter_bcf:
    input:
        expand("results/hard_filtering/subsets/{{bcfsubsamp}}/bisnps-{sgc}-maf-0.05-{{bcfsubsamp}}.bcf", sgc=sg_or_chrom)
    output:
        bcf="results/bcf/{bcfsubsamp}-bisnps-maf-0.05.bcf",
        tbi="results/bcf/{bcfsubsamp}-bisnps-maf-0.05.bcf.csi"
    log:
        "results/logs/hard_filtering/subsets/{bcfsubsamp}/concat-bisnps-maf-0.05-{bcfsubsamp}.log"
    benchmark:
        "results/benchmarks/hard_filtering/subsets/{bcfsubsamp}/concat-bisnps-maf-0.05-{bcfsubsamp}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2>{log}; "

## this rule filters the biallelic snp and maf file to include only autosomal regions and not scaffold reads
rule aut_subset_bcf:
    input:
        bcf="results/bcf/{bcfsubsamp}-bisnps-maf-0.05.bcf",
        tbi="results/bcf/{bcfsubsamp}-bisnps-maf-0.05.bcf.csi",
        regions="results/psmc/get-aut-regions/autosomal_regions.bed"
    output:
        bcf="results/bcf/{bcfsubsamp}-aut-bisnps-maf-0.05.bcf",
        tbi="results/bcf/{bcfsubsamp}-aut-bisnps-maf-0.05.bcf.csi"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/subsets/{bcfsubsamp}/aut-bisnp-maf-0.05-{bcfsubsamp}.log"
    benchmark:
        "results/benchmarks/hard_filtering/subsets/{bcfsubsamp}/aut-bisnp-maf-0.05-{bcfsubsamp}.bmk"
    shell:
        " (bcftools view -Ob -R {input.regions} {input.bcf} > {output.bcf}; "
        " bcftools index {output.bcf}) 2> {log} "

###############################
### end of subset bcf section ###
###############################



#### all callable sites section ####

#### the following is essentially a copy of all the rules above ####
#### but run on the vcf sections including all callable sites, not just variant sites ####

#### GATK section ####
# selects only snps by chrom or scaffold group from the missing-corrected vcf sections
rule make_snp_vcf_all_callable_sites:
    input:
        vcf="results/calling/corrected_missing_vcf_sect/all_callable_sites/{sg_or_chrom}.vcf.gz",
        tbi="results/calling/corrected_missing_vcf_sect/all_callable_sites/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/all_callable_sites/select-variants/snps-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/all_callable_sites/select-variants/snps-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/hard_filtering/all_callable_sites/selectvariants/select-snps-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/all_callable_sites/selectvariants/select-snps-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type SNP -O {output.vcf} > {log} 2>&1 "

#select only indels by chrom or scaffold group from the missing-corrected vcf sections
rule make_indel_vcf_all_callable_sites:
    input:
        vcf="results/calling/corrected_missing_vcf_sect/all_callable_sites/{sg_or_chrom}.vcf.gz",
        tbi="results/calling/corrected_missing_vcf_sect/all_callable_sites/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/all_callable_sites/select-variants/indels-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/all_callable_sites/select-variants/indels-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/hard_filtering/all_callable_sites/selectvariants/select-indels-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/all_callable_sites/selectvariants/select-indels-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type INDEL -O {output.vcf} > {log} 2>&1 "

# filter the snps and indels according to the recommendations at:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
rule hard_filter_snps_all_callable_sites:
    input:
        vcf="results/hard_filtering/all_callable_sites/select-variants/snps-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/all_callable_sites/select-variants/snps-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/all_callable_sites/variant-filtered/snps-filtered-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/all_callable_sites/variant-filtered/snps-filtered-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/hard_filtering/all_callable_sites/variantfiltration/hard-filter-snps-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/all_callable_sites/variantfiltration/hard-filter-snps-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    params:
        filters=config["params"]["gatk"]["VariantFiltration"]["snps"]
    shell:
        "gatk VariantFiltration "
        " -V {input.vcf} "
        " {params.filters} "
        " -O {output.vcf} > {log} 2>&1 "

rule hard_filter_indels_all_callable_sites:
    input:
        vcf="results/hard_filtering/all_callable_sites/select-variants/indels-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/all_callable_sites/select-variants/indels-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/all_callable_sites/variant-filtered/indels-filtered-{sg_or_chrom}.vcf.gz",
        idx="results/hard_filtering/all_callable_sites/variant-filtered/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/hard_filtering/all_callable_sites/variantfiltration/hard-filter-indels-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/all_callable_sites/variantfiltration/hard-filter-indels-{sg_or_chrom}.bmk"
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
rule merge_filtered_vcfs_all_callable_sites:
    input:
        snp="results/hard_filtering/all_callable_sites/variant-filtered/snps-filtered-{sg_or_chrom}.vcf.gz",
        indel="results/hard_filtering/all_callable_sites/variant-filtered/indels-filtered-{sg_or_chrom}.vcf.gz",
        snp_idx="results/hard_filtering/all_callable_sites/variant-filtered/snps-filtered-{sg_or_chrom}.vcf.gz.tbi",
        indel_idx="results/hard_filtering/all_callable_sites/variant-filtered/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/hard_filtering/all_callable_sites/merged-filtered/hardfiltered-merged-{sg_or_chrom}.vcf.gz",
    log:
        "results/logs/hard_filtering/all_callable_sites/merge_filtered_vcfs/hardfiltered-merged-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/all_callable_sites/merge_filtered_vcfs/hardfiltered-merged-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk.yaml"
    shell:
        " gatk MergeVcfs -I {input.indel} -I {input.snp} -O {output.vcf} 2> {log} "

## convert each merged filtered vcf into a bcf
# I think there's a way to do this all in one with concat but I was running into errors
rule merged_filtered_vcfs_to_bcfs_all_callable_sites:
    input:
        vcf="results/hard_filtering/all_callable_sites/merged-filtered/hardfiltered-merged-{sg_or_chrom}.vcf.gz",
    output:
        bcf="results/hard_filtering/all_callable_sites/merged-filtered-bcf/hardfiltered-merged-{sg_or_chrom}.bcf",
    log:
        "results/logs/hard_filtering/all_callable_sites/merged_filtered_vcfs_to_bcfs/hardfiltered-merged-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/all_callable_sites/merged_filtered_vcfs_to_bcfs/hardfiltered-merged-{sg_or_chrom}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools view -Ob {input} > {output} 2>{log} "

###############################


# this gives us a single merged bcf file without any maf filtering & with variants that failed the hard filtering flagged
rule bcf_concat_hardfilter_merged_sect_all_callable_sites:
    input:
        expand("results/hard_filtering/all_callable_sites/merged-filtered-bcf/hardfiltered-merged-{sgc}.bcf", sgc=sg_or_chrom),
    output:
        bcf="results/bcf/all_callable_sites/all_callable_sites.bcf",
        tbi="results/bcf/all_callable_sites/all_callable_sites.bcf.csi"
    log:
        "results/logs/hard_filtering/all_callable_sites/bcf_concat_hardfilter_merged_sect/bcf_concat_log.txt"
    benchmark:
        "results/benchmarks/hard_filtering/all_callable_sites/bcf_concat_hardfilter_merged_sect/bcf_concat.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2> {log}; "

###############################

### this removes the variants that failed the GATK filters ###
### and also removes any variant that is not a biallelic snp ###
### I am not sure what it will do to non-variant sites? ###

## this will run on each scaffold group or chromosome 
# for each of the maf options specified in the maf_cutoff part of the config.yaml file
# I have it set up to run on the missing-corrected-merged-vcf files from the correct_merged_vcfs rule 
rule filter_bisnp_all_callable_sites:
    input:
        "results/hard_filtering/all_callable_sites/merged-filtered/hardfiltered-merged-{sg_or_chrom}.vcf.gz"
    output:
        "results/hard_filtering/all_callable_sites/filter-bisnp/bisnps-{sg_or_chrom}.bcf"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/all_callable_sites/filter_bisnp/bisnps-{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/hard_filtering/all_callable_sites/filter_bisnp/bisnps-{sg_or_chrom}.bmk"
    shell:
        " bcftools view -Ob -v snps -m 2 -M 2 -i 'FILTER=\"PASS\"' "
        " {input} > {output} 2>{log} "

## this gives merged the sections into a bcf
rule concat_bisnp_bcf_all_callable_sites:
    input:
        expand("results/hard_filtering/all_callable_sites/filter-bisnp/bisnps-{sgc}.bcf", sgc=sg_or_chrom)
    output:
        bcf="results/bcf/all_callable_sites/bisnps_all_callable_sites.bcf",
        tbi="results/bcf/all_callable_sites/bisnps_all_callable_sites.bcf.csi"
    log:
        "results/logs/hard_filtering/all_callable_sites/concat_bisnp_bcf/bisnps.log"
    benchmark:
        "results/benchmarks/hard_filtering/all_callable_sites/concat_bisnp_bcf/bisnps.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2>{log}; "

## this rule filters the biallelic snp and maf file to include only autosomal regions and not scaffold reads
rule get_aut_bcf_all_callable_sites:
    input:
        bcf="results/bcf/all_callable_sites/bisnps_all_callable_sites.bcf",
        tbi="results/bcf/all_callable_sites/bisnps_all_callable_sites.bcf.csi",
        regions="results/psmc/get-aut-regions/autosomal_regions.bed"
    output:
        bcf="results/bcf/all_callable_sites/aut_bisnps_all_callable_sites.bcf",
        tbi="results/bcf/all_callable_sites/aut_bisnps_all_callable_sites.bcf.csi"
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/hard_filtering/all_callable_sites/aut_bisnp_bcf/aut_bisnps_all_callable_sites.log"
    benchmark:
        "results/benchmarks/hard_filtering/all_callable_sites/aut_bisnp_bcf/aut_bisnps_all_callable_sites.bmk"
    shell:
        " (bcftools view -Ob -R {input.regions} {input.bcf} > {output.bcf}; "
        " bcftools index {output.bcf}) 2> {log} "

###############################
### end of all callable sites section ### 
###############################







### The following 2 rules are copied from Eric's post-bcf workflow (bcftools_filter.smk) with edits ###
## https://github.com/eriqande/mega-post-bcf-exploratory-snakeflows
## they generate a BCF file with only biallelic snps that pass a MAF cuttoff of 0.05
# I don't currently use these right now (11/12/24)
rule bcf_filt_scatter:
    input:
        bcf="results/bcf/all.bcf",
        csi="results/bcf/all.bcf.csi",
        sfile="config/sample_subsets/all-fish.txt",
        scat_path="results/scatter_config/scatters_1200000.tsv"
    params:
        bcftools_opts="-v snps -m 2 -M 2 -i 'FILTER=\"PASS\" && MAF >= 0.05'",
        scat="{scatter}"
    output:
        scat_regions=temp("results/pca/scat_regions/{scatter}.scat_regions.tsv"),
        bcf=temp("results/bcf/filt_biallelic_maf_0.05/sections/{scatter}.bcf"),
        stats=temp("results/bcf/filt_biallelic_maf_0.05/sections/{scatter}.bcf_stats.txt"),
        pos=temp("results/bcf/filt_biallelic_maf_0.05/sections/{scatter}.positions.tsv.gz"),
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/bcf_filt_scatter/bcf/filt_biallelic_maf_0.05/sections/{scatter}.log"
    benchmark:
        "results/benchmarks/bcf_filt_scatter/filt_biallelic_maf_0.05/sections/{scatter}.bmk"
    shell:
        " (awk -v scat='{params.scat}' -f workflow/scripts/pca/get_scat_regions.awk {input.scat_path} > {output.scat_regions} && "
        "    bcftools view -Ou -R {output.scat_regions} {input.bcf} | "
        "    bcftools +fill-tags -Ou -- -t all | "
        "    bcftools view -Ob {params.bcftools_opts} > {output.bcf}  && "
        "    bcftools stats --af-tag MAF {output.bcf} > {output.stats} && "
        "    bcftools query -f '%CHROM\\t%POS\\n' {output.bcf} | gzip -c > {output.pos} "
        " ) 2> {log} "

rule bcf_filt_gather:
    input:
        bcfs=expand("results/bcf/filt_biallelic_maf_0.05/sections/{scat}.bcf", scat=unique_scats),
        poses=expand("results/bcf/filt_biallelic_maf_0.05/sections/{scat}.positions.tsv.gz", scat=unique_scats),
        statses=expand("results/bcf/filt_biallelic_maf_0.05/sections/{scat}.bcf_stats.txt", scat=unique_scats),
    output:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
        pos="results/bcf/filt_biallelic_maf_0.05/info/positions.tsv.gz",
        stats="results/bcf/filt_biallelic_maf_0.05/info/bcf_stats.txt",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/bcf/filt_biallelic_maf_0.05/main.log"
    benchmark:
        "results/benchmarks/bcf/filt_biallelic_maf_0.05/main.bmk"
    shell:
        " ( bcftools concat --naive {input.bcfs} | "
        "  bcftools sort -Ob > {output.bcf} && "
        "  bcftools index {output.bcf} && "
        "  cat {input.poses} > {output.pos} && "
        "  plot-vcfstats -m {input.statses} > {output.stats} "
        " ) 2> {log} "