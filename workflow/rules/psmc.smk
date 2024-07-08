## rule to get a consensus fastq sequence file for PSMC
# option -C 50 downgrades mapping quality (by coeff given) for reads containing excessive mismatches
# option -d sets and minimum read depth and -D sets the maximum 
# It is recommended to set -d to a third of the average depth and -D to twice
# this takes the bams from each individual and generates a tmp vcf file then 
# generates a consensus sequence for psmc
# alternatively add a rule that does this from the hard filtered vcf files after
# joint calling them in the calling phase
rule psmc_consensus_sequence:
    input:
        bam="results/angsd_bams/overlap_clipped/{sample}.bam", #should this be a bam with -baq 2 done on it for indel stuff? 
        ref="resources/genome/OmykA.fasta",
    output:
        "results/psmc/bams2psmc/psmc-consensus-sequence/{sample}.fq.gz"
    conda:
        "../envs/sambcftools.yaml"
    resources:
        time="23:59:59"
    log:
        "results/logs/psmc/bams2psmc/psmc-consensus-sequence/{sample}.log"
    benchmark:
        "results/benchmarks/psmc/bams2psmc/psmc-consensus-sequence/{sample}.bmk"
    shell:
        "bcftools mpileup -C50 -f {input.ref} {input.bam} | bcftools call -c - | " 
        "vcfutils.pl vcf2fq -d 6 -D 36 | gzip > {output} 2> {log}"