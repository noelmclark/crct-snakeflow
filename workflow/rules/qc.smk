#rule to calculate avg depth of coverage of the mkdup bam file using samtools
rule get_coverage_depth:
    input:
        "results/angsd_bams/overlap_clipped/{sample}.bam"
    output:
        "results/qc/coverage-depth/{sample}.txt"
    conda:
        "../envs/sambcftools.yaml"
    log:
        "results/logs/qc/get-coverage-depth/{sample}.log"
    benchmark:
        "results/benchmarks/qc/get-coverage-depth/{sample}.bmk"
    shell:
        "samtools depth -a -H {input} -o {output} 2> {log}"

rule calc_avg_depth:
    input:
        "results/qc/coverage-depth/{sample}.txt",
    output:
        temp("results/qc/coverage-depth/avg-count/avg-depth-{sample}.txt"),
    log:
        "results/logs/qc/coverage-depth/avg-count/avg-depth-{sample}.log"
    benchmark:
        "results/benchmarks/qc/coverage-depth/avg-count/avg-depth-{sample}.bmk"
    shell:
        " awk '{{sum+=$3}} END {{ print \"has average depth of =\", sum/NR}}' {input} "
        " > {output} 2> {log}"

rule concat_avg_depth:
    input:
        expand("results/qc/coverage-depth/avg-count/avg-depth-{s}.txt", s=sample_list),
    output:
        "results/qc/coverage-depth/avg-count/avg-depth.txt",
    log:
        "results/logs/qc/coverage-depth/avg-count/avg-depth.log"
    benchmark:
        "results/benchmarks/qc/coverage-depth/avg-count/avg-depth.bmk"
    shell:
        " for i in {input}; do "
        " (echo $i && cat $i); done > {output} 2> {log}"

# maybe using this to determine which individual per population to use for hPSMC
rule count_above_10_depth:
    input:
       "results/qc/coverage-depth/{sample}.txt",
    output:
        temp("results/qc/coverage-depth/avg-count/above-10-{sample}.txt"),
    log:
        "results/logs/qc/coverage-depth/avg-count/above-10-{sample}.log",
    benchmark:
        "results/benchmarks/qc/coverage-depth/avg-count/above-10-{sample}.bmk",
    shell:
        " awk '$3 >= 10 {{ count++ }} END {{ print \"has a count of bases with depth greater than 10 =\", count}}' {input}"
        " > {output} 2> {log}"

rule concat_counts:
    input:
       expand("results/qc/coverage-depth/avg-count/above-10-{s}.txt", s=sample_list),
    output:
        "results/qc/coverage-depth/avg-count/above-10.txt",
    log:
        "results/logs/qc/coverage-depth/avg-count/above-10.log",
    benchmark:
        "results/benchmarks/qc/coverage-depth/avg-count/above-10.bmk",
    shell:
        " for i in {input}; do "
        " (echo $i && cat $i); done > {output} 2> {log}"

############

rule samtools_stats:
    input:
        "results/angsd_bams/overlap_clipped/{sample}.bam"
    output:
        "results/qc/samtools_stats/{sample}.txt",
    log:
        "results/logs/qc/samtools_stats/{sample}.log",
    benchmark:
        "results/benchmarks/qc/samtools_stats/{sample}.bmk",
    conda:
        "../envs/sambcftools.yaml"
    shell:
        "samtools stats {input} > {output} 2> {log} "


############
## Rule to calculate % missingness of individual BAM files using a custom bash script and bedtools
rule calc_missingness:
    input:
        "results/psmc/aut-bams/aut_{sample}.bam",
    output:
        temp("results/qc/missingness/{sample}.txt"),
    conda:
        "../envs/bedtools.yaml",
    resources:
        time="6:00:00",
    log:
        "results/logs/qc/missingness/{sample}.log",
    benchmark:
        "results/benchmarks/qc/missingness/{sample}.bmk",
    shell:
        " ./workflow/scripts/PSMC/calc_percent_missingness.sh {input} {output} > {log} "

rule concat_missingness:
    input:
        expand("results/qc/missingness/{s}.txt", s=sample_list),
    output:
        "results/qc/missingness/percent-missing-aut-bams.txt",
    log:
        "results/logs/qc/missingness/percent-missing-aut-bams.log",
    benchmark:
        "results/benchmarks/qc/missingness/percent-missing-aut-bams.bmk",
    shell:
        " for i in {input}; do "
        " (echo $i && cat $i); done > {output} 2> {log}"
############

## get bcftools stats

rule bcftools_stats_autbisnpno5indel:
    input:
        bcf="results/bcf/aut-bisnps-no5indel.bcf",
        csi="results/bcf/aut-bisnps-no5indel.bcf.csi"
    output:
        "results/qc/bcftools_stats/{sample}-stats-aut-bisnps-no5indel.txt",
    log:
        "results/logs/qc/bcftools_stats/{sample}-stats-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/qc/bcftools_stats/{sample}-stats-aut-bisnps-no5indel.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools view -Ou -s {wildcards.sample} {input.bcf} | "
        " bcftools stats > {output} 2> {log} "

rule bcftools_stats_allcallablesites:
    input:
        bcf="results/bcf/all_callable_sites/all_callable_sites.bcf",
        csi="results/bcf/all_callable_sites/all_callable_sites.bcf.csi"
    output:
        "results/qc/bcftools_stats/{sample}-stats-all-callable-sites.txt",
    log:
        "results/logs/qc/bcftools_stats/{sample}-stats-all-callable-sites.log",
    benchmark:
        "results/benchmarks/qc/bcftools_stats/{sample}-stats-all-callable-sites.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools view -Ou -s {wildcards.sample} {input.bcf} | "
        " bcftools stats > {output} 2> {log} "


###########################################################################################
## next 2 copied from Eric
# here is a rule to use bcftools stats to summarize what is found in the
# final BCF files.  Basically, we want to run bcftools stats with the option
# to compile statistics about every single sample, and also to make a histogram
# of our NMISS INFO tag that we made earlier.  The wildcarding here is so that
# we can do this three different ways:  1) with all the variants 2) with only those
# variants that PASS the filter, and 3) with only those variants that don't pass
# the hard filtering. Note that fitering in the bcftools stats command seems to
# conflict with samples in there so that things don't work, so I am piping it
# THIS is done on each chromo or scaff group since doing it on the full file can
# take upwards of 5 hours on a big data set.
rule bcf_section_summaries:
    input:
        "results/hard_filtering/both-filtered-{sg_or_chrom}.bcf",
    output:
        "results/qc/bcftools_stats/sections/{sg_or_chrom}-{filter_condition}.txt",
    log:
        "results/logs/qc/bcftools_stats/{sg_or_chrom}-{filter_condition}.log",
    params:
        comma_samples=",".join(unique_sample_ids),
        filter_opt=get_bcftools_stats_filter_option,
        stop=len(unique_sample_ids),
        steps=len(unique_sample_ids) + 1
    benchmark:
        "results/benchmarks/bcftools_stats/{sg_or_chrom}-{filter_condition}.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view {params.filter_opt} -Ou {input} | "
        " bcftools stats -s {params.comma_samples} "
        " -u NMISS:0:{params.stop}:{params.steps} > "
        " {output} 2> {log} "



# this is messed up.  Too much duplication and bullshit, but I
# want to do the maf-filtered ones, too.
rule bcf_maf_section_summaries:
    input:
        "results/hard_filtering/both-filtered-{sg_or_chrom}-maf-{maf}.bcf",
    output:
        "results/qc/bcftools_stats/maf_sections/{sg_or_chrom}-maf-{maf}.txt",
    log:
        "results/logs/qc/bcftools_stats/{sg_or_chrom}-maf-{maf}.log",
    params:
        comma_samples=",".join(unique_sample_ids),
        stop=len(unique_sample_ids),
        steps=len(unique_sample_ids) + 1
    benchmark:
        "results/benchmarks/qc/bcftools_stats/{sg_or_chrom}-maf-{maf}.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools stats -s {params.comma_samples} "
        " -u NMISS:0:{params.stop}:{params.steps}  {input} > "
        " {output} 2> {log} "


rule combine_bcftools_stats:
    input:
        expand("results/qc/bcftools_stats/sections/{sgc}-{{filter_condition}}.txt", sgc=unique_chromosomes + unique_scaff_groups)
    output:
        "results/qc/bcftools_stats/all-{filter_condition}.txt",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/qc/combine_bcftools_stats/{filter_condition}.log",
    shell:
        " plot-vcfstats -m {input} > {output} 2>{log} "


rule combine_maf_bcftools_stats:
    input:
        expand("results/qc/bcftools_stats/maf_sections/{sgc}-maf-{{maf}}.txt", sgc=unique_chromosomes + unique_scaff_groups)
    output:
        "results/qc/bcftools_stats/all-pass-maf-{maf}.txt",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/qc/combine_maf_bcftools_stats/maf-{maf}.log",
    shell:
        " plot-vcfstats -m {input} > {output} 2>{log} "