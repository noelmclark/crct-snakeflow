#!/usr/bin/env Rscript

.libPaths("/projects/nomclark@colostate.edu/software/anaconda/envs/rCNV/lib/R/library") #set library path to the conda env where I installed R and rCNV
library("rCNV")
library("readr")

# Define the directory containing VCF files
input_dir <- "./results/rCNV-by-scat/vcf"           # Change this to your actual VCF folder path
output_dir <- "./results/rCNV-by-scat/out"  # Change this if you want a different output path
dvs_dir <- "./results/rCNV-by-scat/dvs" 
cnv_dir <- "./results/rCNV-by-scat/cnv" 

# Create output directory if it doesn't exist
if (!dir.exists(dvs_dir)) {
  dir.create(dvs_dir)
}
if (!dir.exists(cnv_dir)) {
  dir.create(cnv_dir)
}

# List all .vcf files in the directory
vcf_files <- list.files(input_dir, pattern = "\\.vcf$", full.names = TRUE)
ainfo_files <- list.files(output_dir, pattern = "\\.tsv$", full.names = TRUE)

# Process each VCF file to generate allele info table
for (vcf_file in vcf_files) { #replace vcf_files with leftover if there are files already
  cat("Processing:", vcf_file, "\n")
  
  # Construct output file name
  base_name <- tools::file_path_sans_ext(basename(vcf_file))
  dvsfile <- file.path(dvs_dir, paste0(base_name, "_dvs.tsv"))
  cnvfile <- file.path(cnv_dir, paste0(base_name, "_cnv.tsv"))
  
  # Load VCF
  vcf <- readVCF(vcf_file)

  # Calculate Fis
  fis <- mean(h.zygosity(vcf)$Fis)
  
  # Load allele info
  allele_file <- file.path(output_dir, paste0(base_name, "_allele_info.tsv"))
  AINFO <- read_tsv(allele_file)

  # Get deviant table
  dvs <- dupGet(AINFO, test = c("z.05","chi.05"), Fis = fis)
  
  # Get CNV table using kmeans clustering method
  CV <- cnv(AINFO, test=c("z.05","chi.05"), filter = "kmeans")
  
  # Save output
  write_tsv(dvs, file = dvsfile)
  write_tsv(CV, file = cnvfile)

  # Remove intermediate files from global environment before restarting
  rm(vcf, fis, AINFO, dvs, CV)
  
  cat("Finished processing and removed from environment:", vcf_file, "\n")
}

cat("All VCF files processed.\n")


## now process dvs and cnv tables to remove flagged variants and join together into one tsv

# Define new directories
dvs_clean_dir <- "./results/rCNV-by-scat/dvs-clean" 
cnv_clean_dir <- "./results/rCNV-by-scat/cnv-clean" 

# Create output directory if it doesn't exist
if (!dir.exists(dvs_clean_dir)) {
  dir.create(dvs_clean_dir)
}
if (!dir.exists(cnv_clean_dir)) {
  dir.create(cnv_clean_dir)
}

dvs_files <- list.files(dvs_dir, pattern = "\\.tsv$", full.names = TRUE)
cnv_files <- list.files(cnv_dir, pattern = "\\.tsv$", full.names = TRUE)

for (dvs_file in dvs_files) {
    cat("Cleaning:", dvs_file, "\n")

    # Construct output file name
    base_name <- tools::file_path_sans_ext(basename(dvs_file))
    dvs_clean <- file.path(dvs_clean_dir, paste0(base_name, "_dvs_clean.tsv"))
    
    # Load the table
    dvs <- read_tsv(dvs_files)
    
    # Filter dvs files by removing rows with "deviant" flag
    dvs_cleaned <- dvs %>% filter(dup.stat == "non-deviant")
    
    # write new file 
    write_tsv(dvs_cleaned, file = dvs_clean)

    # Remove intermediate files from global environment before restarting
    rm(dvs, dvs_cleaned)

    cat("Finished cleaning and saved:", dvs_file, "\n")
}