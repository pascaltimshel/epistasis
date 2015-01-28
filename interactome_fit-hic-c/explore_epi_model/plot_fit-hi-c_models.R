############### SYNOPSIS ###################
# CLASS: script; understanding; investigation of pvalues; visualization

# PURPOSE: *PLOT* *fit-hi-c* models.

# Genotypes: EGCUT clean (832 individuals, ~2e6 SNPs)
# Phenotype/probe: ~10,000 top variance/mean probes (egcut.peer_residuals_log2_k50.top50_mean_top50_var_refseq.txt)

############################################

library(ggplot2)
library(reshape2)
library(snpStats)

############
rm(list=ls())
wd <- "~/git/epistasis/interactome_fit-hic-c/explore_epi_model"
setwd(wd)
############


######################################## Read probes ###########################################
### Read file with 9269 probes and 832 individuals
#file.probes <- "/Users/pascaltimshel/Dropbox/Ubuntu/peer_analysis/export_residuals/egcut.peer_residuals_log2_k50.top50_mean_top50_var_refseq.txt"
# df.probes <- read.table(file.probes, h=T)
# rownames(df.probes) <- df.probes$IID
# df.probes <- subset(df.probes, select=-c(IID, FID))
# save(df.probes, file="RData_tmp/tmp_EGCUT_9269_probes.RData")
load(file="/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/R_epistasis_tests/RData_tmp/tmp_EGCUT_9269_probes.RData") # df.probes



######################################## Read Genotypes ###########################################
#### Small test set
bed.file <- "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/SNPs12959/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.fit-hi-c.SNPs12959"
#### Read data:
sample <- read.plink(bed.file)
# Converting SnpMatrix to Numeric Matrix
geno <- as(sample$genotypes, "numeric") # 0,1,2. Missing data will be NA.

######################################## File checks ###########################################
df.check.sync.probe_geno <- data.frame(rownames(df.probes), rownames(geno))
all(df.check.sync.probe_geno[,1]==df.check.sync.probe_geno[,2]) # --> TRUE: in sync

######################################## Read FastEpistasis results ###########################################
# HEADER: CHR  SNP_A	CHR.1	SNP_B	BETA	CHISQ	PVALUE	PHENOTYPE
cols.fastEpi <- c("chr1", "snp1", "chr2", "snp2", "BETA", "CHISQ", "PVALUE", "probename")
file.fastEpi <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/fastEpi_1e-10_0108_212539/fastEpi_table.p_value.sorted.top1000.csv"
df.fastEpistasis.results <- read.csv(file=file.fastEpi, col.names=cols.fastEpi, stringsAsFactors=F)
str(df.fastEpistasis.results)

### Create ID column: probe.snp1.snp2
df.fastEpistasis.results$ID <- apply(df.fastEpistasis.results[,c("probename", "snp1", "snp2")], 1, function(df.row) {paste(df.row[1], paste(sort(c(as.character(df.row[2]), as.character(df.row[3]))),collapse="."),sep=".")})
str(df.fastEpistasis.results)


##################################### Plot epistasis model using ###########################################
source("function_plot_epistasis_model.R")
idx <- 1:10
#idx <- 1:nrow(df.snp_probe_pairs.subset)
plots <- plot_EpiModel(idx, df.fastEpistasis.results, add_data_heatmap=TRUE, save_images=TRUE, save_significant_treshold=1e-5)

df.fastEpistasis.results[1, "snp1"]

