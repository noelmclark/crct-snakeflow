#!/usr/bin/env Rscript

.libPaths("/projects/nomclark@colostate.edu/software/anaconda/envs/rCNV/lib/R/library") #set library path to the conda env where I installed R and rCNV
library("rCNV")
library("readr")

#load the paths from the snakemake rule
vcf_file <- snakemake@input$vcf
ainfo_file <- snakemake@input$ainfo
dvs_file <- snakemake@output$dvs
cnv_file <- snakemake@output$cnv

# Load VCF
vcf <- readVCF(vcf_file)

# Calculate Fis
fis <- mean(h.zygosity(vcf)$Fis)

# Load allele info tsv
AINFO <- read_tsv(ainfo_file)

# Get deviant table
dvs <- dupGet(AINFO, test = c("z.05","chi.05"), Fis = fis)

# Get CNV table using kmeans clustering method
CV <- cnv(AINFO, test=c("z.05","chi.05"), filter = "kmeans")

# Save output
write_tsv(dvs, file = dvs_file)
write_tsv(CV, file = cnv_file)
