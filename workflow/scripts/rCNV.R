#redirect output and messages/erros to the log
#log <- file(snakemake@log[[1]], open="wt")
#sink(log, type = "output")
#sink(log, type = "message")
if (!require("rCNV")) install.packages("rCNV", repos = "http://cran.us.r-project.org")

library(rCNV)

#load the paths from the snakemake rule
vcf_file <- snakemake@input$vcf
outfile <- snakemake@output$tsv

vcf<-readVCF(vcf_file)

fis<-mean(h.zygosity(vcf)$Fis)
gt<-hetTgen(vcf,"GT")
ad<-hetTgen(vcf,"AD")
ad[ad==""] <- "0,0"
out<-allele.info(X=ad, x.norm=gt, Fis=fis)
# X is the corrected non-normalized allele depth table and x.norm is the normalized allele depth table
# out is a data frame with duplication status of each snp, and other stats.

#write the output to the file path
write_tsv(out, file = outfile)