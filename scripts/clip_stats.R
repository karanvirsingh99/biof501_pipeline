library(data.table)
library(tidyverse)

setwd("data/summary_stats/")

AN <- fread(paste0("../../", snakemake@input[[1]]))
BMI <- fread(paste0("../../", snakemake@input[[2]]))
thousandgenomes <- read_table("../reference_genomes/reference.1000G.maf.0.005.txt.gz")
chr_three <- subset(thousandgenomes, CHR==3)

print("Files loaded")

AN_ch3 <- subset(AN, ID %in% chr_three$SNP)
AN_ch3 <- AN_ch3 %>% rename(A1 = "ALT",
                            A2 = "REF")
BMI_chr3 <- subset(BMI, SNP %in% chr_three$SNP)

write_tsv(AN_ch3, file=gzfile("AN_subset.tsv.gz"))
write_tsv(BMI_chr3, file=gzfile("BMI_subset.tsv.gz"))
