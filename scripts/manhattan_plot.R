##################################################################
##                        Manhattan Plot                        ##
##################################################################
library(qqman)
library(tidyverse)

# Clean up data
load(snakemake@input[[1]])
non_BMI_AN <- outputGWAS[[2]]
non_BMI_AN <- non_BMI_AN %>% rename(P = "Pval_Estimate")

pdf("results/non_BMI_AN_manhattan_CHR3.pdf", width=12, height=10)
non_BMI_AN_manhattan <- manhattan(non_BMI_AN, annotatePval = 5E-8, annotateTop = FALSE)
dev.off()
