library(data.table)
library(GenomicSEM)

# STEP 1: Clean up GWAS summary statistics

# - AN summary statistics from Watson et al.
AN <- fread(file=snakemake@input[[1]])

# - BMI summary statistics from
BMI <- fread(file=snakemake@input[[2]])

#################################################################
##               Calculate effective sample size               ##
#################################################################

#========
##  AN   =
##========

#read in information about sample size per cohort
AN_cohort<-read.csv("data/summary_stats/an_gwas_cohorts.csv")

#calculate sample prevalence for each cohort
AN_cohort$v<-AN_cohort$Post.QC_Ncases/(AN_cohort$Post.QC_Ncases+AN_cohort$Post.QC_Ncontrols)

#calculate cohort specific effective sample size
AN_cohort$EffN<-4*AN_cohort$v*(1-AN_cohort$v)*(AN_cohort$Post.QC_Ncases+AN_cohort$Post.QC_Ncontrols)

#calculate sum of effective sample size: 46321.9
effective_sample_size_AN <- sum(AN_cohort$EffN)

# This is unneccessary since the *NEFFDIV2* column in the summary stats is just the effective sum sample size divided by 2!!!

#========
##  BMI  =
##========

BMI_sample_size <- 322154

#################################################################
##                         munge stats                         ##
#################################################################
setwd("./results")
AN_file <- paste0("../", snakemake@input[[1]])
BMI_file <- paste0("../", snakemake@input[[2]])

munge(AN_file,
      "../data/reference_genomes/w_hm3.noMHC.snplist",
      trait.names="AN",
      N=effective_sample_size_AN,
      info.filter = 0.7,
      maf.filter = 0.01)

munge(BMI_file,
      "../data/reference_genomes/w_hm3.noMHC.snplist",
      trait.names="BMI",
      N=322154,
      info.filter = 0.7,
      maf.filter = 0.01)


#################################################################
##                     LD score regression                     ##
#################################################################

traits <- c("AN.sumstats.gz","BMI.sumstats.gz")
sample.prev <- c(NA,NA)
population.prev <- c(NA,NA)
ld<-"../data/reference_genomes/eur_w_ld_chr/"
wld <- "../data/reference/genomes/eur_w_ld_chr/"
trait.names<-c("AN", "BMI")

LDSCoutput <- ldsc(traits,
                   sample.prev,
                   population.prev,
                   ld,
                   wld,
                   trait.names)


save(LDSCoutput, file="LDSCoutput.RData")

##################################################################
##                    SEM model without SNPs                    ##
##################################################################

load(file="LDSCoutput.RData")

model<-'B=~NA*AN + start(0.4)*BMI
        NB=~NA*AN

         NB~~1*NB
         B~~1*B
         B~~0*NB

         BMI ~~ 0*AN
         BMI~~0*BMI
         AN~~0*AN'


output<-usermodel(LDSCoutput,estimation="DWLS",model=model)


