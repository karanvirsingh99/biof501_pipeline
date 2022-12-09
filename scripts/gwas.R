##################################################################
##                        Model with SNP                        ##
##################################################################
library(GenomicSEM)
setwd("results")


load(paste0("../", snakemake@input[[1]]))
load(paste0("../", snakemake@input[[2]]))

# model with SNP
model<-'B=~NA*AN +start(0.2)*AN + start(0.5)*BMI
         NB=~NA*AN +start(0.2)*AN

         B~SNP
         NB~SNP

         NB~~1*NB
         B~~1*B
         B~~0*NB

         BMI ~~ 0*AN
         BMI~~0*BMI
         AN~~0*AN
         SNP~~SNP'

#Run the Genomic SEM GWAS
outputGWAS<-userGWAS(covstruc=LDSCoutput,
                     SNPs=p_sumstats,
                     estimation="DWLS",
                     model=model,sub =c("B~SNP","NB~SNP"))# printwarn = FALSE

save(outputGWAS, file="outputGWAS.RData")
load("outputGWAS.RData")
