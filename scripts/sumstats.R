library(GenomicSEM)
setwd("results")

#################################################################
##               How to choose sumstat arguments               ##
#################################################################

# ------ AN ------

# Case/control study using logistic regression.
# The standard errors are for an odds ratio
# The BETA column is the ln(OR), or a logistic beta

# set se.logit=T, OLS=F, linprob=F

# ------ BMI ------

# Continuous trait
# set se.logit=F, OLS=T, linprob= F


files = c(paste0("../", snakemake@input[[1]]), paste0("../", snakemake@input[[2]]))
ref = "../data/reference_genomes/reference.1000G.maf.0.005.txt.gz"
trait.names = c("AN","BMI")
se.logit = c(T,F)
info.filter = 0.7
maf.filter = 0.01

p_sumstats<-sumstats(files, ref,trait.names, se.logit, info.filter, maf.filter, OLS=c(F,T),
                                 linprob=NULL, N=c(46321.90,322154), betas = NULL)

save(p_sumstats, file="Sumstats.RData")



