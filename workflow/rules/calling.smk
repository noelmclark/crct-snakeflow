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

