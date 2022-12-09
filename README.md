# Snakemake workflow: GWAS-by-subtraction

### By: Karanvir Singh (Dennis Lab)

------------------------------------------------------------------------

# Project description

This workflow fits a genomic structural equation model (GenomicSEM) using genome-wide association study (GWAS) summary statistics for anorexia nervosa ([Watson et. al](https://www.nature.com/articles/s41588-019-0439-2)) and BMI ([Locke et. al](https://www.nature.com/articles/nature14177)). The model regresses the effect of each variant identified in the GWASs on two latent variables - BMI and non-BMI. This multivariate effectively computes the effect size of genetic variants identified that are associated with anorexia independently of BMI.

## Overview

Anorexia nervosa (AN) has the highest mortality out of all psychiatric conditions, occurs early in life and is often chronic (1). No effective medication for AN exists, and only approximately 50% of AN patients achieve full recovery (2). Individuals who experience the psychological symptoms of AN but do not present with a subthreshold body weight (which is often the case in children) are diagnosed with atypical AN (AAN). Prognosis of AN and AAN can be improved by earlier detection and treatment. However, we know very little about the antecedents of AN, and especially AAN, which hampers early intervention. Genes account for 48-74% of AN risk (3), while the genetic risk for AAN has not yet been studied. A recent genome wide association study (GWAS) identified variants that additively explain between 11-17% of AN phenotypic variation. (4) In addition, a major finding of the GWAS of AAN was that some of the genetic variants were also associated with metabolic traits (4). This paints a framework in which the genetic risk of AN can be separated between a metabolic and psychological component (mirroring the clinical diagnosis of the disorder). This finding sets the stage for leveraging the existing GWAS of AN in parallel with existing GWAS of other metabolic traits to identify the genetics of AAN. In this pipeline, I use genomic structural equation modeling to model the genetic risk for AN that is independent of BMI. This method has been previously used to model the genetic architecture of educational achievement that is independent from the genetics of cognitive potential. (5)

The model can be visualized like this:

![](images/GWAS%20by%20substraction-02.png)

There are three observed variables: a SNP, BMI and AN, and two latent variables: BMI and non-BMI. The model estimates two paths for each SNP: one that corresponds to the SNP's effect on AN through BMI, while the second corresponds to the SNP's effect on AN that is not correlated with BMI at all. In this pipeline, we estimate the latter for all SNPs on chromosome 3 to identify significant variants that affect AN genetic risk independently of BMI.

# Workflow Overview

This workflow is based on the gwas-by-subtraction tutorial available [here](https://rpubs.com/MichelNivard/565885)

Here is a visualization of the workflow:

![](images/Screen%20Shot%202022-12-09%20at%202.49.58%20PM.png)

The key steps of this workflow are:

1.  Munge summary statistics and run linkage disequilibrium score regression

This step filters SNPs in our summary statistics based on minor allele frequency and low quality (INFO). It then computes the genetic covariance matrix between AN and BMI, and the standard errors associated with the entries in this matrix. The ldsc.log file contains useful information about the SNP-based heritability (h2) for both traits, as well as the genetic correlation. The .RData file will be used to build our model later in the pipeline.

2.  Clip summary statistics

To maintain the runtime of this workflow within reasonable limits, I subset the summary statistics for both traits to only include SNPs in the chromosome three before fitting the model.

3.  Sumstats

Analagous to Step 1, this step prepares the summary statistics by aligning them to a reference genomes (1000 Genomes Phase 3), merge them together, and transform the effect sizes before fitting the model. This step will output a AN_BMI_sumstats.log file, which documents how many SNPs are present in each file before and after merging and how columns are interpreted by the software. The Sumstats.RData output from this step is used to fit the model in Step 4.

4.  Fit the gwas-by-subtraction model

This step fits the gwas-by-subtraction model, using the sumstats from step 3 and the LDSC output from step 1. There are two outputs from this step: a outputGWAS.RData file, which can be imported to explore results further. This file is a list that contains our two GWAS outputs: the non-BMI AN GWAS (effect sizes on AN that are independent of BMI), and a BMI-AN GWAS (effect sizes on AN that depend on BMI). The non-BMI AN GWAS results are exported as non_BMI_AN.txt

5.  Plot manhattan plot

The final step of this pipeline generates a manhattan plot for chromosome 3 based on the non-BMI-AN GWAS result. Any significant hits at the genome wide level (p\< 5e-8) are highlighted.

## Dependencies

Follow instructions to install [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) , [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation) and [mamba](https://mamba.readthedocs.io/en/latest/installation.html) before running this pipeline.

This snakemake workflow relies on the following main package dependencies:

``` sh
snakemake
python==3.6
r-base==4.2.2
r-tidyverse
r-data.table
r-qqman
r-genomicsem
```

## Installation

Deactivate any conda environments you might have running by opening the terminal and running this command

``` sh
conda deactivate
```

To access the workflow, clone this repository to your computer by running:

``` sh
git clone git@github.com:karanvirsingh99/biof501_pipeline.git
```

Navigate to the project directory (change the project directory to where you cloned the repository)

``` sh
cd biof501_pipeline
```

The repository contains an environment.yml file which mamba can use to install all dependencies the pipeline needs to run. You can run the following command to build this environment, which will be named biof501_pipeline.

``` sh
mamba env create -f environment.yml
```

Activate the environment:

``` sh
conda activate biof501_pipeline
```

### Setup data

Some of the files needed for the pipeline are too large to be stored on github.

From the project directory, run the following command:

``` sh
bash get_files.sh
```

This will grab the LD-score reference, AN GWAS summary statistics and place them in the correct folders.

Additionally, the 1000 Genomes Phase 3 reference file needs to be downloaded from [here](https://drive.google.com/file/d/10MbM3JZW44hthxqX7gaVRr_UU75Z05DS/view?usp=sharing). Place this file in data/reference_genomes

## Usage

Once all dependencies are installed, and you have activated the environment, you can run the pipeline by running this command, as long as you are in the project directory where the Snakefile is located. Please note that Step 4 is the most computationally intensive and should take approximately thirty minutes to run. Refer to the expected outputs to ensure that all steps were completed successfully.

``` sh
snakemake --cores 1
```

## Expected outputs

*Step 1 - munge and ldsc*

| File                                    | Intepretation                                                                                  |
|------------------------------------|------------------------------------|
| AN_munge.log                            | Logs how columns in the summary statistics file for AN were interpreted by the munge function  |
| BMI_munge.log                           | Logs how columns in the summary statistics file for BMI were interpreted by the munge function |
| AN.sumstats.gz_BMI.sumstats.gz_ldsc.log | LDSC summary results for AN and BMI, such as SNP-based heritability and genetic correlation    |
| LDSCoutput.RData                        | Covariance matrix that will be used in Step 4                                                  |

*Step 2 - subsetting summary statistics*

| File              | Intepretation                                      |
|-------------------|----------------------------------------------------|
| AN_subset_tsv.gz  | AN summary statistics with only chromosome 3 SNPs  |
| BMI_subset_tsv.gz | BMI summary statistics with only chromosome 3 SNPs |

*Step 3 - sumstats*

| File                | Intepretation                                                                                                                                                                                  |
|------------------------------------|------------------------------------|
| AN_BMI_sumstats.log | Log of sumstats function - useful to check how columns were interpreted, how many SNPs are lost due to MAF, INFO filter, overlap between reference genome, and overlap between the two traits. |
| Sumstats.RData      | Harmonized summary statistics with transformed effect sizes, which will be used in Step 4                                                                                                      |

*Step 4 - gwas-by-subtraction model fitting*

| File             | Intepretation                                                                                                                         |
|------------------------------------|------------------------------------|
| outputGWAS.RData | List with two entries: the first contains the BMI-dependent effect sizes on AN, and the second the BMI-independent-effect sizes on AN |

*Step 5 - manhattan plot*

| File                          | Intepretation                                                                                                                                                               |
|------------------------------------|------------------------------------|
| non_BMI_AN_manhattan_CHR3.pdf | A manhattan plot for chromosome 3, plotting the BMI-independent effect size on AN for all chr3 SNPs. Variants significant at the genome-wide level (p\<5e-8) are annotated, |

## Manhattan plot

This is the final expected output for the pipeline, a manhattan plot for all SNPs in chromosome 3, with p-values based on their effect size on AN that is independent of BMI. All SNPs significant at the genome-wide level (p \< 5E-8) are annotated.

![](images/Screen%20Shot%202022-12-09%20at%2012.40.39%20PM.png)

## Next steps

Results from this pipeline can be further investigated. Significant hits identified in chromosome 3 can be fine-mapped using tools like [plink](https://www.cog-genomics.org/plink/2.0/) to identify the most likely causal variants. Functional mapping using colocalization with known QTLs can also be used to determine the effect of SNPs in non-coding regions, as well as to identify associated target genes.

	<h1>Directory Tree</h1><p>
	This is what your directory should look like after running the pipeline
	<a href=".">.</a><br>
	├── <a href="./README.md">README.md</a><br>
	├── <a href="./Snakefile">Snakefile</a><br>
	├── <a href="./data/">data</a><br>
	│   ├── <a href="./data/reference_genomes/">reference_genomes</a><br>
	│   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/">eur_w_ld_chr</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/1.l2.M_5_50">1.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/1.l2.ldscore.gz">1.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/10.l2.M_5_50">10.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/10.l2.ldscore.gz">10.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/11.l2.M_5_50">11.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/11.l2.ldscore.gz">11.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/12.l2.M_5_50">12.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/12.l2.ldscore.gz">12.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/13.l2.M_5_50">13.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/13.l2.ldscore.gz">13.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/14.l2.M_5_50">14.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/14.l2.ldscore.gz">14.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/15.l2.M_5_50">15.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/15.l2.ldscore.gz">15.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/16.l2.M_5_50">16.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/16.l2.ldscore.gz">16.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/17.l2.M_5_50">17.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/17.l2.ldscore.gz">17.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/18.l2.M_5_50">18.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/18.l2.ldscore.gz">18.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/19.l2.M_5_50">19.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/19.l2.ldscore.gz">19.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/2.l2.M_5_50">2.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/2.l2.ldscore.gz">2.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/20.l2.M_5_50">20.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/20.l2.ldscore.gz">20.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/21.l2.M_5_50">21.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/21.l2.ldscore.gz">21.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/22.l2.M_5_50">22.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/22.l2.ldscore.gz">22.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/3.l2.M_5_50">3.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/3.l2.ldscore.gz">3.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/4.l2.M_5_50">4.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/4.l2.ldscore.gz">4.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/5.l2.M_5_50">5.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/5.l2.ldscore.gz">5.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/6.l2.M_5_50">6.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/6.l2.ldscore.gz">6.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/6_old.l2.M_5_50">6_old.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/6_old.l2.ldscore.gz">6_old.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/7.l2.M_5_50">7.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/7.l2.ldscore.gz">7.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/8.l2.M_5_50">8.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/8.l2.ldscore.gz">8.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/9.l2.M_5_50">9.l2.M_5_50</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/9.l2.ldscore.gz">9.l2.ldscore.gz</a><br>
	│   │   │   ├── <a href="./data/reference_genomes/eur_w_ld_chr/README">README</a><br>
	│   │   │   └── <a href="./data/reference_genomes/eur_w_ld_chr/w_hm3.snplist">w_hm3.snplist</a><br>
	│   │   ├── <a href="./data/reference_genomes/reference.1000G.maf.0.005.txt.gz">reference.1000G.maf.0.005.txt.gz</a><br>
	│   │   └── <a href="./data/reference_genomes/w_hm3.noMHC.snplist">w_hm3.noMHC.snplist</a><br>
	│   └── <a href="./data/summary_stats/">summary_stats</a><br>
	│   &nbsp;&nbsp;&nbsp; ├── <a href="./data/summary_stats/AN_subset.tsv.gz">AN_subset.tsv.gz</a><br>
	│   &nbsp;&nbsp;&nbsp; ├── <a href="./data/summary_stats/BMI_subset.tsv.gz">BMI_subset.tsv.gz</a><br>
	│   &nbsp;&nbsp;&nbsp; ├── <a href="./data/summary_stats/SNP_gwas_mc_merge_nogc.tbl.uniq.gz">SNP_gwas_mc_merge_nogc.tbl.uniq.gz</a><br>
	│   &nbsp;&nbsp;&nbsp; ├── <a href="./data/summary_stats/an_gwas_cohorts.csv">an_gwas_cohorts.csv</a><br>
	│   &nbsp;&nbsp;&nbsp; └── <a href="./data/summary_stats/pgcAN2.2019-07.vcf.tsv.gz">pgcAN2.2019-07.vcf.tsv.gz</a><br>
	├── <a href="./environment.yml">environment.yml</a><br>
	├── <a href="./get_files.sh">get_files.sh</a><br>
	├── <a href="./images/">images</a><br>
	│   ├── <a href="./images/GWAS%20by%20substraction-02.png">GWAS by substraction-02.png</a><br>
	│   ├── <a href="./images/Screen%20Shot%202022-12-09%20at%2012.40.39%20PM.png">Screen Shot 2022-12-09 at 12.40.39 PM.png</a><br>
	│   └── <a href="./images/Screen%20Shot%202022-12-09%20at%202.49.58%20PM.png">Screen Shot 2022-12-09 at 2.49.58 PM.png</a><br>
	├── <a href="./results/">results</a><br>
	│   ├── <a href="./results/AN.sumstats.gz">AN.sumstats.gz</a><br>
	│   ├── <a href="./results/AN.sumstats.gz_BMI.sumstats.gz_ldsc.log">AN.sumstats.gz_BMI.sumstats.gz_ldsc.log</a><br>
	│   ├── <a href="./results/AN_BMI_sumstats.log">AN_BMI_sumstats.log</a><br>
	│   ├── <a href="./results/AN_munge.log">AN_munge.log</a><br>
	│   ├── <a href="./results/BMI.sumstats.gz">BMI.sumstats.gz</a><br>
	│   ├── <a href="./results/BMI_munge.log">BMI_munge.log</a><br>
	│   ├── <a href="./results/LDSCoutput.RData">LDSCoutput.RData</a><br>
	│   ├── <a href="./results/Sumstats.RData">Sumstats.RData</a><br>
	│   ├── <a href="./results/non_BMI_AN_manhattan_CHR3.pdf">non_BMI_AN_manhattan_CHR3.pdf</a><br>
	│   └── <a href="./results/outputGWAS.RData">outputGWAS.RData</a><br>
	└── <a href="./scripts/">scripts</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="./scripts/clip_stats.R">clip_stats.R</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="./scripts/gwas.R">gwas.R</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="./scripts/manhattan_plot.R">manhattan_plot.R</a><br>
	&nbsp;&nbsp;&nbsp; ├── <a href="./scripts/munge_and_ldsc.R">munge_and_ldsc.R</a><br>
	&nbsp;&nbsp;&nbsp; └── <a href="./scripts/sumstats.R">sumstats.R</a><br>
<br><br><p>

## References

1.  van Eeden, A. E., van Hoeken, D. & Hoek, H. W. Incidence, prevalence and mortality of anorexia nervosa and bulimia nervosa. Curr. Opin. Psychiatry **34**, 515--524 (2021).
2.  Steinhausen, H.-C. Outcome of eating disorders. Child Adolesc. Psychiatr. Clin. N. Am. **18**, 225--242 (2009).
3.  Yilmaz, Z., Hardaway, J. A. & Bulik, C. M. Genetics and Epigenetics of Eating Disorders. Adv. Genomics Genet. **5**, 131--150 (2015).
4.  Watson, H. J. et al. Genome-wide association study identifies eight risk loci and implicates metabo-psychiatric origins for anorexia nervosa. Nat. Genet. **51**, 1207--1214 (2019)
5.  Demange, P. A. *et al.* Investigating the genetic architecture of noncognitive skills using GWAS-by-subtraction. *Nat Genet* **53**, 35--44 (2021).

