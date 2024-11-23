## this rule counts the total variable sites needed before counting hets in next rule
rule bcftools_stats:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf",
        csi="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf.csi",
    output:
        "results/inbreeding/het/{sample}-stats.txt",
    conda:
        "../envs/bcftools.yaml",
    log:
        "results/logs/inbreeding/het/{sample}-stats.log",
    benchmark:
        "results/benchmarks/inbreeding/het/{sample}-stats.bmk"
    shell:
        " bcftools view -Ou -s {wildcards.sample} {input.bcf} | "
        " bcftools stats > {output} 2> {log} "

## This rule counts the heterozygous sites (0/1) in an individual vcf
rule count_hets_nmiss:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf",
        csi="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf.csi",
    output:
        het="results/inbreeding/het/{sample}-het-count.txt",
        nmiss="results/inbreeding/het/{sample}-nmiss-count.txt",
    conda:
        "../envs/bcftools.yaml",
    log:
        "results/logs/inbreeding/het/{sample}-het-nmiss-count.log",
    benchmark:
        "results/benchmarks/inbreeding/het/{sample}-het-nmiss-count.bmk"
    shell:
        " ( bcftools view -H -s {wildcards.sample} {input.bcf} | grep \"0/1:\" | wc -l > {output.het} && "
        " bcftools view -H -s {wildcards.sample} {input.bcf} | grep -v \"./.:\" | wc -l > {output.nmiss} ) "
        " 2> {log} "

rule calc_percent_het:
    input:
        het="results/inbreeding/het/{sample}-het-count.txt",
        nmiss="results/inbreeding/het/{sample}-nmiss-count.txt",
    output:
        "results/inbreeding/het/{sample}-het-perc.txt",
    log:
        "results/logs/inbreeding/het/{sample}-het-perc.log",
    benchmark:
        "results/benchmarks/inbreeding/het/{sample}-het-perc.bmk"
    shell:
        " awk 'NR==FNR{{het=$0; next}} {{print het / $0}}' {input.het} {input.nmiss} "
        " > {output} 2> {log} "


# thank you to chat gpt for this code
rule combine_het_perc:
    input:
        results=expand("results/inbreeding/het/{s}-het-perc.txt", s=sample_list),
    output:
        "results/inbreeding/het/combined-het-percents.tsv",
    log:
        "results/logs/inbreeding/het/combined-het-percents.tsv"
    benchmark:
        "results/benchmarks/inbreeding/het/combined-het-percents.tsv"
    shell:
        """
        echo -e "Sample\tPercHet" > {output}
        for file in {input.results}; do
            sample=$(basename $file -het-perc.txt)
            value=$(cat $file)
            echo -e "$sample\t$value" >> {output}
        done
        """

##################################


## another way to get het from plink gt counts using awk
rule get_het_from_gt_count:
    input:
        gcount="results/plink/gt-count/{sample}-aut-bisnps-no5indel.gcount",
    output:
        "results/inbreeding/het/from-plink-gt-counts/{sample}-aut-bisnps-no5indel.txt",
    log:
        "results/logs/inbreeding/het/from-plink-gt-counts/{sample}-aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/inbreeding/het/from-plink-gt-counts/{sample}-aut-bisnps-no5indel.bmk",
    shell:
        """
        awk 'BEGIN {{FS="\\t"; OFS="\\t"; sum=0; count=0}} 
        NR>1 && $10 == 0 {{sum+=$6; count++}} 
         END {{
             if (count > 0) print "Sum of HET_REF_ALT_CTS:", sum; 
             print "Total variants:", count; 
             print "Average HET_REF_ALT_CTS:", (count ? sum/count : 0)
         }}' {input.gcount} > {output} 2> {log}
        """

## this rule takes the outputs from above and combines them into a tsv
rule get_het_from_gt_count:
    input:
        dir="results/inbreeding/het/from-plink-gt-counts/",
    output:
        "results/inbreeding/het/aut-bisnps-no5indel.tsv",
    conda:
        "../envs/hpsmc-split.yaml" #bc it already has python in it
    log:
        "results/logs/inbreeding/het/aut-bisnps-no5indel.log",
    benchmark:
        "results/benchmarks/inbreeding/het/aut-bisnps-no5indel.bmk",
    shell:
        " python het_from_gt_counts.py {input.dir} {output} "




































## Not currently used

# this filters the bcf file by sample names listed in the given population wildcard 
# and generates an allele frequency file for each population
rule make_bcftools_pop_afreq:
    input:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
    params:
        pops=get_comma_sep_pop_names
    output:
        afreq="results/roh/allele-freq/snps-maf-0.05/{population}-freqs.tabs.gz",
        tbi="results/roh/allele-freq/snps-maf-0.05/{population}-freqs.tabs.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    log:
        freq="results/logs/roh/allele-freq/snps-maf-0.05/{population}-freqs.log",
        index="results/logs/roh/allele-freq/snps-maf-0.05/{population}-freqs.tbi.log"
    benchmark:
        "results/benchmarks/roh/allele-freq/snps-maf-0.05/{population}-freqs.bmk",
    shell:
        " bcftools view -s {params.pops} {input.bcf} |"
        " bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' | "
        " bgzip -c > {output.afreq} 2> {log.freq} && "
        " tabix -f -s1 -b2 -e2 {output.afreq} > {output.tbi} 2> {log.index} "


rule run_bcftools_roh:
    input:
        bcf="results/bcf/filt_biallelic_maf_0.05/main.bcf",
        csi="results/bcf/filt_biallelic_maf_0.05/main.bcf.csi",
        afreq="results/roh/allele-freq/snps-maf-0.05/{population}-freqs.tabs.gz",
        tbi="results/roh/allele-freq/snps-maf-0.05/{population}-freqs.tabs.gz.tbi"
    params:
        pops=get_comma_sep_pop_names
    output:
        "results/roh/snps-maf-0.05/{population}-roh.txt",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/roh/snps-maf-0.05/{population}-roh.log",
    benchmark:
        "results/benchmarks/roh/snps-maf-0.05/{population}-roh.bmk",
    shell:
        " bcftools view -s {params.pops} {input.bcf} | "
        " bcftools roh --AF-file {input.afreq} --GTs-only 30 "
        " -o {output} 2> {log} "