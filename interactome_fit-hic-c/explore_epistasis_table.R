############### SYNOPSIS ###################
# CLASS: script; analysis; explorative
# PURPOSE: read and explore epistasis_table

# *IMPORTANT*: this script originated as a MODIFIED DUPLICATE of another script.
# VERSION v.1 --> "compare_R_and_FastEpistasis_p-values_inflated_p_values.R" from the "..XXX/EGCUT_DATA/R_epistasis_tests" directory
# VERSION v.2 --> "analyze_FastEpistasis_p-values_fit-hi-c.R"

# STEPS
# 1) Read *epistasis_table*
# 2a) Extract SNPs and write reduced PLINK files
# 2b) Read PLINK files (this is needed for MAF calculation and R-anova model)
# 3) Calculate MAF
# 4) Run RunEpistasisTests
# 5) Make plots
# 6) OPTIONAL: export merged results

############################################

library(snpStats)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tools)

rm(list=ls())

wd <- path.expand("~/git/epistasis/interactome_fit-hic-c")
setwd(wd)
######################


################################################################################################
###################################### Settings/params ########################################
################################################################################################

### Job
#param.job_identifier <- "hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8" # n=38,987
param.job_identifier <- "hESC_width_2500_maf_5_q_1e-13_epi1_1e-10" # n=248,510

### Epistasis table
param.epistasis_table.type <- "epistasis_table_pruned_hemani.txt" # this is needed for loading the correct epistasis table

### Paths - out
path.out <- path.expand("~/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2")
path.out.base <- file.path(path.out, param.job_identifier)
path.out.base.analysis_epi_table <- file.path(path.out.base, "analysis_epi_table")
path.out.base.analysis_epi_table

path.out.base.plots_model <- file.path(path.out.base.analysis_epi_table, "plot_epistasis_model")

path.out.base.PLINK_files <- path.out.base.analysis_epi_table

### Paths - libary (scripts)
path.library <- path.expand("~/git/epistasis/interactome_fit-hic-c/explore_epi_model")

### Files
file.epi_table <- file.path(path.out.base, param.epistasis_table.type)
file.epi_table # e.g. "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8/epistasis_table_pruned_hemani.txt"



################################################################################################
######################################## Initialize ##########################################
################################################################################################

### Source ###
source(file.path(path.library,"function_write_snp_subset_plink_file.R"))
  # write_snp_subset_plink_file()
source(file.path(path.library,"function_subset_min_genotype_class_count.R"))
  # subset_genotype_class_count()
source(file.path(path.library,"function_Replication_fastEpistasis_comparison.R"))
  # RunEpistasisTests()
source(file.path(path.library,"function_plot_epistasis_model.R"))
  # plot_EpiModel()


### Create directories if needed
if (file.exists(path.out.base.analysis_epi_table)) {
  print("OBS: path.out.base.analysis_epi_table folder exists")
}
dir.create(path.out.base.analysis_epi_table, showWarnings=FALSE) # dir.create() does NOT crash if the directory already exists, it just prints out a warning.
print(path.out.base.analysis_epi_table) # no trailing slash


################################# LOAD epistasis_table ############################
### EXAMPLE - April 2015 ###
### REMARKS
  # - no header
# 5  rs7722032  21	rs7281390	-0.23753	57.45355	3.46056E-14	ILMN_1699160	hic_1_5965	True
# 1	rs11122320	2	rs4954710	+0.12256	53.48095	2.61113E-13	ILMN_1688234	hic_1_2276	True
# ...
# 4  rs7676504	18	rs1943313	+0.33703	58.93330	1.63108E-14	ILMN_1791759	null_1_7165	True
# 12	rs7966152	11	rs968064	+0.32100	53.49671	2.59026E-13	ILMN_2065022	null_1_2091	True

### READ table
#epi_table.header <- c("CHR","SNP_A","CHR","SNP_B","BETA","CHISQ","PVALUE","PHENOTYPE","EIID","SIGNIFICANCE") # KEEP THIS HEADER UPDATED!
epi_table.header <- c("chr1", "snp1", "chr2", "snp2", "BETA", "CHISQ", "PVALUE", "probename","EIID","bonferroni_significant") # KEEP THIS HEADER UPDATED!
df.epi_table <- read.delim(file.epi_table, stringsAsFactors=F, header=F, col.names=epi_table.header)#
str(df.epi_table)
nrow(df.epi_table)

### ***Creating EID column***
df.epi_table$EID <- sub("_\\d*?$", "", df.epi_table$EIID, perl=TRUE) # OBS: meta characters need to be escaped (e.g. \d --> \\d)
# REF for R regex: http://www.johndcook.com/blog/r_language_regex/

### Overview of EID
df.overview <- df.epi_table %>% group_by(EID) %>% summarise(count=n()) %>% arrange(desc(count))
df.overview

### UNIQUE SNPs
SNP_AB.union <- union(df.epi_table$snp1, df.epi_table$snp2); length(SNP_AB.union)



######################################## Read probes ###########################################
### Read file with 9269 probes and 832 individuals
#file.probes <- "/Users/pascaltimshel/Dropbox/Ubuntu/peer_analysis/export_residuals/egcut.peer_residuals_log2_k50.top50_mean_top50_var_refseq.txt"
# df.probes <- read.table(file.probes, h=T)
# rownames(df.probes) <- df.probes$IID
# df.probes <- subset(df.probes, select=-c(IID, FID))
# save(df.probes, file="RData_tmp/tmp_EGCUT_9269_probes.RData")
load(file="/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/R_epistasis_tests/RData_tmp/tmp_EGCUT_9269_probes.RData") # df.probes

###############################################################################################################
######################################## Read FastEpistasis results ###########################################
###############################################################################################################

### *BACKWARDS COMPATIBILITY* ###
# We are using the data frame with variable name "df.fastEpistasis.results" because of legacy of previously written code
df.fastEpistasis.results <- df.epi_table

### Create ID column: probe.snp1.snp2
df.fastEpistasis.results$ID <- apply(df.fastEpistasis.results[,c("probename", "snp1", "snp2")], 1, function(df.row) {paste(df.row[1], paste(sort(c(as.character(df.row[2]), as.character(df.row[3]))),collapse="."),sep=".")})
str(df.fastEpistasis.results)

###############################################################################################################
######################################## Write PLINK files ###########################################
###############################################################################################################

### WRITING NEW SET OF PLINK FILES (from function_write_snp_subset_plink_file.R)
# --> output set of plink files will be written to PATH.OUT/egcut_snpsubset.{bed,bim,fam,log}
df.snpstats <- write_snp_subset_plink_file(df=df.fastEpistasis.results, path.out=path.out.base.PLINK_files)
df.snpstats

######################################## Read Genotypes ###########################################
#### Small test set
#bed.file <- "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/SNPs12959/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.fit-hi-c.SNPs12959"
#bed.file <- "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/all_clean/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean"
bed.file <- paste0(path.out.base.PLINK_files, "/egcut_snpsubset")
#### Read data:
sample <- read.plink(bed.file)
# Converting SnpMatrix to Numeric Matrix
geno <- as(sample$genotypes, "numeric") # 0,1,2. Missing data will be NA.
dim(geno)

######################################## File checks ###########################################
df.check.sync.probe_geno <- data.frame(rownames(df.probes), rownames(geno))
all(df.check.sync.probe_geno[,1]==df.check.sync.probe_geno[,2]) # --> TRUE: in sync



#################################### Calculate min class count ######################################

source(file.path(path.library,"function_subset_min_genotype_class_count.R"))
# From script: function_subset_min_genotype_class_count.R
#idx <- 1:1000
idx <- 1:nrow(df.fastEpistasis.results)
df.output <- subset_genotype_class_count(idx=idx, df.input=df.fastEpistasis.results, min_genotype_class_count_threshold=3)
nrow(df.output)
  # hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8 | min 3 --> n=703

# [1] "plotting #71266/#248510"
# [1] "plotting #71267/#248510"
# Show Traceback
# 
# Rerun with Debug
# Error in data.frame(snp1.allele, snp2.allele) : 
#   arguments imply differing number of rows: 832, 0 

df.overview.2 <- df.output %>% group_by(EID) %>% summarise(count=n()) %>% arrange(desc(count))
df.overview.2

######################################## Calculate MAF ###########################################
### snpStats - MAF for all SNPs in the loaded genotype file
snpsum <- col.summary(sample$genotypes)
df.snpsum.maf <- data.frame(snp=rownames(snpsum), MAF=snpsum$MAF)
head(df.snpsum.maf)

### Adding column to data frame: MAF
df.fastEpistasis.results$snp1.maf <- df.snpsum.maf$MAF[match(df.fastEpistasis.results$snp1, df.snpsum.maf$snp)]
df.fastEpistasis.results$snp2.maf <- df.snpsum.maf$MAF[match(df.fastEpistasis.results$snp2, df.snpsum.maf$snp)]
### Alternative approach: LOOP OVER VECTOR (sapply()) - MUCH SLOWER!
#snp1.maf <- sapply(df.fastEpistasis.results$snp1, function(rsID) {df.snpsum.maf$MAF[match(rsID, df.snpsum.maf$snp)]})
str(df.fastEpistasis.results)

### Adding column to data frame: pairwise MAF
df.fastEpistasis.results$snp.maf.min <- with(df.fastEpistasis.results, pmin(snp1.maf, snp2.maf))
#sum(df.fastEpistasis.results$snp.maf.min > 0.2)

#df.fastEpistasis.results.maf.gt <- subset(df.fastEpistasis.results, snp.maf.min>=0.04)
df.fastEpistasis.results.maf.gt <- subset(df.fastEpistasis.results, snp.maf.min>=0.2)
df.fastEpistasis.results.maf.gt

######################################## Run RunEpistasisTests ###########################################
df.probe_snp_pair <- subset(df.fastEpistasis.results.maf.gt, select=c(snp1,snp2, probename))
#df.probe_snp_pair <- subset(df.fastEpistasis.results, select=c(snp1,snp2, probename))
#df.probe_snp_pair <- subset(df.fastEpistasis.results[1:100,], select=c(snp1,snp2, probename))
str(df.probe_snp_pair)

### From script "function_Replication_fastEpistasis_comparison.R"
# df.probe_snp_pair: 3 columns. snp1;snp2,probename
df.res <- RunEpistasisTests(idx=1:nrow(df.probe_snp_pair), 
                            df.probe_snp_pair=df.probe_snp_pair, 
                            geno=geno, 
                            expression=df.probes)
str(df.res)

### Create ID column: probe.snp1.snp2
df.res$ID <- apply(df.res[,c("probename", "snp1", "snp2")], 1, function(df.row) {paste(df.row[1], paste(sort(c(as.character(df.row[2]), as.character(df.row[3]))),collapse="."),sep=".")})
## OLD METHOD USING sapply (2x slower than apply)
#df.res$ID <- sapply(1:nrow(df.res), function(i) {paste(df.res$probename[i], paste(sort(c(as.character(df.res$snp1[i]), as.character(df.res$snp2[i]))),collapse="."),sep=".")})

### MERGE data frames - keep all observations in "df.res" (?)
df.res.merge <- merge(df.res, df.fastEpistasis.results, by=c("ID"), all.x=T)
#df.res.merge <- merge(df.res, df.fastEpistasis.results, by=c("snp1", "snp2"), all.x=T) #  by.x=c("snp1", "snp2"), by.y=c("SNP_A", "SNP_B")


################# Make some plots ####################
ggplot(df.res.merge, aes(x=beta, y=BETA)) + geom_point() + geom_abline() #beta_plot-6x6
#ggsave("beta_plot-6x6.pdf")
ggplot(df.res.merge, aes(x=F_statistic, y=CHISQ)) + geom_point() # F-statistic_plot-6x6
ggplot(df.res.merge, aes(x=-log10(F_pval), y=-log10(PVALUE))) + geom_point() # P-value_plot-6x6
ggplot(df.res.merge, aes(x=-log10(F_pval), y=-log10(PVALUE))) + geom_point(aes(color=n_usable_data_points)) + geom_abline() # P-value_plot_with_abline_and_color-6x6
#ggsave("P-value_plot_with_abline_and_color-6x6.pdf")

### *NEW* ####
## Scatter plot; color=MAF
df.res.merge.edit <- df.res.merge
df.res.merge.edit$min_maf <- pmin(df.res.merge.edit$snp1_maf, df.res.merge.edit$snp2_maf) # pmin()
ggplot(df.res.merge.edit, aes(x=-log10(F_pval), y=-log10(PVALUE))) + geom_point(aes(color=min_maf)) + geom_abline() # 
#ggsave("P-value_plot_with_abline_and_color-minMAF6x6.pdf")

## Histplot MAF # ****THIS IS HERE I STOPPED ****** #
#ggplot(df.res.merge.edit, aes(x=min_maf)) + geom_histogram() + labs(title="MIN MAF, N=1000")
#ggsave("min_maf_n1000.pdf")


################# Export ####################
df.res.merge.export <- df.res.merge
### OLD WAY (before fit-hi-c) probename.x
#df.res.merge.export <- subset(df.res.merge, select=-c(ID.x, probename.y, ID.y))
### CURRENT WAY (using merge 'by=ID')
#write.csv(df.res.merge.export, file="df.res.merge.export.csv")



##################################### Plot epistasis model  ###########################################
idx <- 1:10
#idx <- 1:nrow(df.snp_probe_pairs.subset)
### Run function from script "function_plot_epistasis_model.R"
plots <- plot_EpiModel(idx, df.probe_snp_pair, path.out=path.out.base.plots_model, add_data_heatmap=TRUE, save_images=TRUE, save_significant_treshold=1e-5)
## TODO plot_EpiModel()
# 1) add MAF(SNP1), MAF(SNP2) and min.maf to plot



######################################## Run individual ###########################################

str(df.fastEpistasis.results)
df.fastEpistasis.results[1,]
str(geno)

i = 1
name.probe <- "ILMN_2065022"
name.snp1 <- "rs17421337" 
name.snp2 <- "rs7237309" 


# Extract data
snp1 <- geno[, colnames(geno) == name.snp1]
snp2 <- geno[, colnames(geno) == name.snp2]
probe <- df.probes[, colnames(df.probes) == name.probe]
#cbind(snp1, snp2)
#probe <- rnorm(nrow(geno)) # random phenotype


### MAF calculation
# REMEMBER to remove NA when using sum()
# number_of_individuals_with_snp/number_of_non_missing_datapoints
# number_of_SNP[x]_alleles/total_number_of_alleles
tmp.snp1.maf <- (sum(snp1==1,na.rm=T)*1 + sum(snp1==2, na.rm=T)*2)/(2*sum(!is.na(snp1)))
tmp.snp2.maf <- (sum(snp2==1,na.rm=T)*1 + sum(snp2==2, na.rm=T)*2)/(2*sum(!is.na(snp2)))
snp1.maf <- ifelse(tmp.snp1.maf <= 0.5, tmp.snp1.maf, 1-tmp.snp1.maf)
snp2.maf <- ifelse(tmp.snp2.maf <= 0.5, tmp.snp2.maf, 1-tmp.snp2.maf)
snp1.maf; snp2.maf

tab <- table(data.frame(snp1, snp2))
tab
snp1.allele <- factor(snp1,levels=c("0", "1", "2", NA), exclude=NULL)
snp2.allele <- factor(snp2,levels=c("0", "1", "2", NA), exclude=NULL)
tab.allele <- table(data.frame(snp1.allele, snp2.allele))
tab.allele
df.twolocus <- as.data.frame(tab.allele)
rownames(df.twolocus) <- paste("count_SNP1/SNP2=", df.twolocus$snp1.allele, "/", df.twolocus$snp2.allele, sep="")
df.twolocus <- subset(df.twolocus, select=c(-snp1.allele, -snp2.allele))
df.twolocus.t <- t(df.twolocus)
df.twolocus.t

cor(snp1, snp2, use="pairwise.complete.obs")
length(tab) # replication_nclass
min(tab, na.rm=T) #replication_minclass
sum(!is.na(snp1) & !is.na(snp2)) # replication_nid

## Full
fullmod.multiplicative <- lm(probe ~ snp1 + snp2 + snp1:snp2)
fullmod.multiplicative
summary(fullmod.multiplicative)
## Marginal
margmod.multiplicative <- lm(probe ~ snp1 + snp2)
summary(margmod.multiplicative)
## Test for model reduction
inttest_multiplicative <- anova(margmod.multiplicative, fullmod.multiplicative)
inttest_multiplicative




