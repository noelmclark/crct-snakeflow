#### These rules use Plink2.0 & Plink1.9 on our filtered BCF ####

### PLINK based rules ###
## these rules take our filtered BCF and run PLINK2 on it to generate different outputs
## --allow-extra-chromosomes lets the bed file include non-human chroms (PLINK is defaulted to humans)
## I removed the Y chrom, mitogenome, and unassembled scaffolds using bcftools, but you could use the --not-chr flag in PLINK

## this rule creates population level allele freq count files which I then combine in a custom script
# the --pheno option splits our bcf into one per population that includes all indivdiuals with the population id from the popfile
# this is so --freq can generate population specific allele frequencies vs ones based on the all sample bcf
rule calc_pop_allele_freq:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf.csi",
        popfile="config/plink-popfile.tsv",
    output:
        afreq="results/plink/allele-freq/pop-counts/aut-snps-{maf}",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/allele-freq/pop-counts/aut-snps-{maf}.log",
    benchmark:
        "results/benchmarksplink/allele-freq/pop-counts/aut-snps-{maf}.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --pheno {input.popfile} "
        " --allow-extra-chr "
        " --freq counts --loop-cats population "
        " --out {output.afreq} 2> {log} "

# this rule keeps all the populations together and generates a whole dataset allele freq file
# this is necessary for running PCA with all samples 
rule calc_allele_freq:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf.csi",
    output:
        afreq="results/plink/allele-freq/aut-snps-{maf}",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/allele-freq/aut-snps-{maf}.log",
    benchmark:
        "results/benchmarks/plink/allele-freq/aut-snps-{maf}.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --freq "
        " --out {output.afreq} 2> {log} "

# this rule creates a list of variants that pass the ld prune filter (20kb window, slide by 1kb, r^2 0.4)
# and a --geno 0.1 missingness filter 
# --bad-ld is needed to force plink to run ld pruning with less than 50 samples (we have 48)
rule prune_linkage_disequilibrium:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf.csi",
    output:
        ld="results/plink/ld-prune/aut-snps-{maf}-ld-pruned",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/ld-prune/aut-snps-{maf}-ld-pruned.log",
    benchmark:
        "results/benchmarks/plink/ld-prune/aut-snps-{maf}-ld-pruned.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 " 
        " --indep-pairwise 20kb 0.4 "
        " --bad-ld "
        " --out {output.ld} 2> {log} "

## This rules generates a PCA using Plink2.0 from our filtered BCF
# the --geno 0.01 applies a 10% missingness threshold filter 
rule make_plink_pca:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf.csi",
        afreq="results/plink/allele-freq/aut-snps-0.05.afreq"
    output:
        pca="results/plink/pca/aut-snps-{maf}-pca",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/pca/aut-snps-{maf}-pca.log",
    benchmark:
        "results/benchmarks/plink/pca/aut-snps-{maf}-pca.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --pca "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --read-freq {input.afreq} "
        " --out {output.pca} 2> {log} "


## This rules generates a PCA using Plink2.0 from our filtered BCF
# and using a LD pruned variant set 
# the --geno 0.01 applies a 10% missingness threshold filter that should be redundant when using the new purned sites 
# the --make-bed file also produced the input needed for running ADMIXTURE
rule make_plink_pruned_pca:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf.csi",
        afreq="results/plink/allele-freq/aut-snps-0.05.afreq",
        ld="results/plink/ld-prune/aut-snps-{maf}-ld-pruned.prune.in"
    output:
        pca="results/plink/pca/aut-snps-{maf}-pruned-pca",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/pca/aut-snps-{maf}-pruned-pca.log",
    benchmark:
        "results/benchmarks/plink/pca/aut-snps-{maf}-pruned-pca.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --extract {input.ld} "
        " --geno 0.1 "
        " --read-freq {input.afreq} "
        " --make-bed "
        " --pca "
        " --out {output.pca} 2> {log} "



## This rule calculate pairwise Fst values using the Weir & Cockerham (1984) method 
# on our hard filtered BCF file with only biallelic snps that pass a MAF cutoff of 0.05
# the --geno 0.01 applies a 10% missingness threshold filter that should be redundant when using the new purned sites 
rule make_pw_fst_snp:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf",
        tbi="results/bcf/autosomal-biallelic-snps-maf-{maf}.bcf.csi",
        popfile="config/plink-popfile.tsv",
        ld="results/plink/ld-prune/aut-snps-{maf}-ld-pruned.prune.in"
    output:
        fst="results/plink/pw-fst/aut-snps-{maf}-fst",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/pw-fst/aut-snps-{maf}-fst.log",
    benchmark:
        "results/benchmarks/plink/pw-fst/aut-snps-{maf}-fst.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --const-fid 0 "
        " --geno 0.1 "
        " --extract {input.ld}"
        " --within {input.popfile} population"
        " --fst population method=wc "
        " --out {output.fst} 2> {log} "



## This rule generates a Phylip file from the filtered BCF which is used as input for the IQ Tree and splits-tree programs
# the --geno 0.01 applies a 10% missingness threshold filter that should be redundant when using the new purned sites 
# need to double check file options and recommended filters in PLINK and IQTree
rule make_phylip:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf",
        csi="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf.csi",
        popfile="config/plink-popfile.tsv",
        ld="results/plink/ld-prune/aut-snps-{maf}-ld-pruned.prune.in"
    output:
        phylip="results/plink/phylip/aut-snps-{maf}",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/phylip/aut-snps-{maf}.log",
    benchmark:
        "results/benchmarks/plink/phylip/aut-snps-{maf}.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --extract {input.ld} "
        " --snps-only "
        " --export phylip used-sites "
        " --out {output.phylip} 2> {log} "



## the next two rules calculate genome wide heterozygosity on the full variant set and the ld-pruned one
# --het computes observed and expected homozygous/heterozygous genotype counts for each sample, 
# and reports method-of-moments F coefficient estimates (i.e. (1 - (<observed het. count> / <expected het. count>)))
rule calc_het:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf",
        csi="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf.csi",
        afreq="results/plink/allele-freq/aut-snps-0.05.afreq",
    output:
        het="results/plink/het/aut-snps-{maf}",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/het/aut-snps-{maf}.log",
    benchmark:
        "results/benchmarks/plink/het/aut-snps-{maf}.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --read-freq {input.afreq} "
        " --het "
        " --out {output.het} 2> {log} "

rule calc_pruned_het:
    input:
        bcf="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf",
        csi="results/bcf/autosomal-biallelic-snps-maf-0.05.bcf.csi",
        afreq="results/plink/allele-freq/aut-snps-0.05.afreq",
        ld="results/plink/ld-prune/aut-snps-{maf}-ld-pruned.prune.in"
    output:
        het="results/plink/het/aut-snps-{maf}-pruned",
    conda:
        "../envs/plink.yaml"
    log:
        "results/logs/plink/het/aut-snps-{maf}-pruned.log",
    benchmark:
        "results/benchmarks/plink/het/aut-snps-{maf}-pruned.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --allow-extra-chr "
        " --geno 0.1 "
        " --read-freq {input.afreq} "
        " --extract {input.ld} "
        " --het "
        " --out {output.het} 2> {log} "



### BONEYARD (i.e. rules not used currently) ###

## This rule identifies runs of homozygosity (ROH) using the PLINK1.9 method
# https://www.cog-genomics.org/plink/1.9/ibd#homozyg

## This rule calculate pairwise Fst values using the Weir & Cockerham (1984) method 
# on our hard filtered BCF file with all variants
# populations are loaded as categorical phenotypes first from a popfile
# if we want to do jacknife blocksize as well: https://groups.google.com/g/plink2-users/c/vr8rHzYVhZo/m/syHgeWLYAAAJ
# I don't call for this currently 
rule make_pw_fst_all:
    input:
        bcf="results/bcf/all.bcf",
        tbi="results/bcf/all.bcf.csi",
        popfile="config/plink-popfile.tsv",
    output:
        fst="results/plink/pw-fst/fst-all",
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=18800,
        cpus=2,
    threads: 2
    log:
        "results/logs/plink/pw-fst/fst-all.log",
    benchmark:
        "results/benchmarks/plink/pw-fst/fst-all.bmk",
    shell:
        " plink2 --bcf {input.bcf} "
        " --set-missing-var-ids @:#[b37]\$r,\$a "
        " --new-id-max-allele-len 387 "
        " --allow-extra-chr "
        " --const-fid 0 --within {input.popfile} population"
        " --out {output.fst} --fst population method=wc 2> {log} "