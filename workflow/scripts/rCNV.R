#!/usr/bin/env Rscript

# Define the directory containing VCF files
input_dir <- "./results/rCNV-by-scat/vcf"           # Change this to your actual VCF folder path
output_dir <- "./results/rCNV-by-scat/"  # Change this if you want a different output path

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List all .vcf files in the directory
vcf_files <- list.files(input_dir, pattern = "\\.vcf$", full.names = TRUE)

# Process each VCF file
for (vcf_file in vcf_files) {
  cat("Processing:", vcf_file, "\n")
  
  # Construct output file name
  base_name <- tools::file_path_sans_ext(basename(vcf_file))
  outfile <- file.path(output_dir, paste0(base_name, "_allele_info.tsv"))
  
  # Load VCF
  vcf <- readVCF(vcf_file)

  # Calculate Fis and other stats
  fis <- mean(h.zygosity(vcf)$Fis)
  gt <- hetTgen(vcf, "GT")
  ad <- hetTgen(vcf, "AD")
  ad[ad == ""] <- "0,0"
  adnor <- cpm.normal(ad, method = "MedR")
  
  # Generate allele info
  out <- allele.info(X = ad, x.norm = adnor, Fis = fis)

  # Save output
  write_tsv(out, file = outfile)
}

cat("All VCF files processed.\n")
