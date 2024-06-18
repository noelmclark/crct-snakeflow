# copied from Eric's calling.smk

## the following 2 rules are used to create the chromo and scaffold group interval lists
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
