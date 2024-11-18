## rules to run admixture

rule admixture_first_k:
    input:
        "results/plink/pca/aut-snps-0.05-pruned-pca.bed"
    output:
        "results/admixture/aut-snps-0.05-pruned-pca.10.Q",
    conda:
        "../envs/admixture.yaml"
    log:
        "results/logs/admixture/aut-snps-0.05-pruned-pca.10.log"
    benchmark:
        "results/benchmarks/admixture/aut-snps-0.05-pruned-pca.10.bmk"
    shell:
        "admixture --cv {input} 10"
