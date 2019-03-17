############### SYNOPSIS ###################
# CLASS: script; thesis illustration
# PURPOSE: analyze probe co-localization

# *IMPORTANT*: this script is a PARTIAL COPY of another script "explore_epistasis_table.R"

############################################

library(snpStats)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr) # REMEMBER 
library(tools)

rm(list=ls())

wd <- path.expand("~/git/epistasis/thesis_plots")
setwd(wd)
######################


################################################################################################
###################################### Settings/params ########################################
################################################################################################

### Job
param.job_identifier <- "hIMR90_width_1000_maf_5_q_1e-06_epi1_1e-10" # n=96,083 | EIID_probe
#param.job_identifier <- "hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8" # n=102,339 | EIID_probe

### Epistasis table
# this is needed for loading the correct epistasis table
#param.epistasis_table.type <- "epistasis_table.txt"
#param.epistasis_table.type <- "epistasis_table_pruned_hemani.txt"
param.epistasis_table.type <- "epistasis_table_pruned_EIID_probe.txt"

### Paths - out
path.out.subdir <- "epistasis_table_processing"

path.out <- path.expand("~/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2")
path.out.base <- file.path(path.out, param.job_identifier, path.out.subdir)
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
source(file.path(path.library,"function_min_genotype_class_count.R"))
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
# 5  rs7722032  21  rs7281390	-0.23753	57.45355	3.46056E-14	ILMN_1699160	hic_1_5965	True
# 1	rs11122320	2	rs4954710	+0.12256	53.48095	2.61113E-13	ILMN_1688234	hic_1_2276	True
# ...
# 4  rs7676504	18	rs1943313	+0.33703	58.93330	1.63108E-14	ILMN_1791759	null_1_7165	True
# 12	rs7966152	11	rs968064	+0.32100	53.49671	2.59026E-13	ILMN_2065022	null_1_2091	True

### READ table
#epi_table.header <- c("CHR","SNP_A","CHR","SNP_B","BETA","CHISQ","PVALUE","PHENOTYPE","EIID","SIGNIFICANCE") # KEEP THIS HEADER UPDATED!
#epi_table.header <- c("chr1", "snp1", "chr2", "snp2", "BETA", "CHISQ", "PVALUE", "probename","EIID","EID",) # KEEP THIS HEADER UPDATED!
epi_table.header <- c("chr1", "snp1", "chr2", "snp2", "BETA", "CHISQ", "PVALUE", "probename","EIID","EID",
                      "RefSeq_ID","Symbol","Chromosome","Probe_Chr_Orientation","Probe_Coordinates_start","Probe_Coordinates") # KEEP THIS HEADER UPDATED!
df.epi_table <- read.delim(file.epi_table, stringsAsFactors=F, header=F, col.names=epi_table.header)#
str(df.epi_table)
nrow(df.epi_table)

### ***Creating EID column*** 
#[*OUTCOMMENT ME LATER*]
# df.epi_table$EID <- sub("_\\d*?$", "", df.epi_table$EIID, perl=TRUE) # OBS: meta characters need to be escaped (e.g. \d --> \\d)
## REF for R regex: http://www.johndcook.com/blog/r_language_regex/


### UNIQUE SNPs
SNP_AB.union <- union(df.epi_table$snp1, df.epi_table$snp2); length(SNP_AB.union)


######################################## *CREATE df.enrichment* ###########################################
EIDs <- c("hic_1", paste0("null_",seq(1,1000)))
df.enrichment <- data.frame(EID=EIDs)
#df.enrichment
### Calculate step enrichment
df.enrichment.1.no_min_geno <- df.epi_table %>% group_by(EID) %>% dplyr::summarise(count=n()) %>% arrange(desc(count))
df.enrichment.1.no_min_geno <- dplyr::rename(df.enrichment.1.no_min_geno, count1_no_min_geno = count) # rename variable
df.enrichment.1.no_min_geno
### UPDATE df.enrichment
df.enrichment <- full_join(df.enrichment, df.enrichment.1.no_min_geno, by="EID")



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

###############################################################################################
################################ *** PROBE co-localization *** ################################
###############################################################################################

#load(file="df.min_genotype_class_count.full.hIMR90_l-1000.RData") # df.min_genotype_class_count.full
load(file="df.min_genotype_class_count.full.hIMR90_l-2500.RData") # df.min_genotype_class_count.full


### "All" co-localization
x <- with(df.min_genotype_class_count.full, (Chromosome == chr1) | (Chromosome == chr2))
sum(x)
length(x)
sum(x)/length(x)*100

#hIMR90-l=1000: 9482/96083 | 9.87 %
#hIMR90-l=2500: 9902/102339 | 9.68 %

### "Spatial pairs" co-localization
x <- with(df.min_genotype_class_count.full, (EID == "hic_1") & ((Chromosome == chr1) | (Chromosome == chr2)))
sum(x)
y <- sum(with(df.min_genotype_class_count.full, (EID == "hic_1"))); y
sum(x)/y*100

#hIMR90-l=1000: 17/189 | 8.99 %
#hIMR90-l=2500: 11/149 | 7.38 %

### "All" co-localization | AFTER minGCC >= 3
df.min_genotype_class_count.minGCC <- subset(df.min_genotype_class_count.full, min_genotype_class_count >= 3)
x <- with(df.min_genotype_class_count.minGCC, (Chromosome == chr1) | (Chromosome == chr2))
sum(x)
length(x)
sum(x)/length(x)*100

#hIMR90-l=1000: 106/1169 | 9.07 %
#hIMR90-l=2500: 95/1051 | 9.04 %


### "Spatial pairs" co-localization | AFTER minGCC >= 3
df.min_genotype_class_count.minGCC.hic <- subset(df.min_genotype_class_count.full, (min_genotype_class_count >= 3) & (EID=="hic_1"))
x <- with(df.min_genotype_class_count.minGCC.hic, (Chromosome == chr1) | (Chromosome == chr2))
sum(x)
length(x)
sum(x)/length(x)*100

#hIMR90-l=1000: 1/2 | 50.00 %
#hIMR90-l=2500: 1/2 | 50.00 %


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
df.fastEpistasis.results.maf.gt <- subset(df.fastEpistasis.results, snp.maf.min>=0.1)
#df.fastEpistasis.results.maf.gt <- subset(df.fastEpistasis.results, snp.maf.min>=0.25)
nrow(df.fastEpistasis.results.maf.gt)


#################################### Calculate min class count ######################################

source(file.path(path.library,"function_min_genotype_class_count.R"))
# From script: function_subset_min_genotype_class_count.R
df.input <- df.fastEpistasis.results.maf.gt
#df.input <- df.fastEpistasis.results
#idx <- 1:1000
idx <- 1:nrow(df.input)
df.min_genotype_class_count.full <- min_genotype_class_count(idx=idx, df.input=df.input)
nrow(df.min_genotype_class_count.full)
# hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8 | min 3 --> n=703

### *THRESHOLD* ###
min_genotype_class_count_threshold <- 3 # 3 | 5 | 10

### SUBSET DATA FRAME
df.min_genotype_class_count.sub <- subset(df.min_genotype_class_count.full, min_genotype_class_count >= min_genotype_class_count_threshold)
print(sprintf("number of SNP-pairs satisfying criteria: %s (%.2f %%)", nrow(df.min_genotype_class_count.sub), nrow(df.min_genotype_class_count.sub)/nrow(df.min_genotype_class_count.full)*100))


