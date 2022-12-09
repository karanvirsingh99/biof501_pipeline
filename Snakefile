rule all:
    input:
        "results/non_BMI_AN_manhattan_CHR3.pdf"

rule clip_stats:
    input:
        "data/summary_stats/pgcAN2.2019-07.vcf.tsv.gz",
	"data/summary_stats/SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
    output:
        "data/summary_stats/AN_subset.tsv.gz",
	"data/summary_stats/BMI_subset.tsv.gz"
    script:
        "scripts/clip_stats.R"

rule munge_and_ldsc:
    input:
        "data/summary_stats/pgcAN2.2019-07.vcf.tsv.gz",
	"data/summary_stats/SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
    output:
        "results/LDSCoutput.RData",
	"results/AN.sumstats.gz",
	"results/BMI.sumstats.gz",
	"results/BMI_munge.log",
	"results/AN_munge.log"
    script:
        "scripts/munge_and_ldsc.R"

rule sumstats:
    input:
        "data/summary_stats/AN_subset.tsv.gz",
	"data/summary_stats/BMI_subset.tsv.gz"
    output:
        "results/AN_BMI_sumstats.log",
	"results/Sumstats.RData"
    script:
        "scripts/sumstats.R"

rule gwas:
    input:
        "results/LDSCoutput.RData",
	"results/Sumstats.RData"
    output:
        "results/outputGWAS.RData"
    script:
        "scripts/gwas.R"

rule manhattan:
    input:
        "results/outputGWAS.RData"
    output:
        "results/non_BMI_AN_manhattan_CHR3.pdf"
    script:
        "scripts/manhattan_plot.R"