############### SYNOPSIS ###################
# CLASS: script; analysis
# PURPOSE: ...

############################################

library(reshape2)
library(ggplot2)
library(dplyr)
library(tools)

rm(list=ls())

wd <- path.expand("~/git/epistasis/simulations")
setwd(wd)
############################################


### HapMapII imputed
file.maf_distribution.imputed <- "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/frq_files/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean.frq.cut5.txt"
df.maf_distribution.imputed <- read.table(file.maf_distribution.imputed, h=T)

downsample <- 3.1e5
df.maf_distribution.imputed.sub <- sample_n(df.maf_distribution.imputed, downsample) # sample specific number
#df.maf_distribution.imputed.sub <- sample_frac(df.maf_distribution.imputed, 0.1) # sample fraction

p <- ggplot(df.maf_distribution.imputed.sub, aes(x=MAF)) + geom_histogram(binwidth=0.01)
p <- p + labs(title=sprintf("EGCUT HapMap II Imputed SNPs. Downsample=%d [N_SNPs]", downsample))
p
ggsave(file="maf_distribution_EGCUT_imputed.pdf" ,w=8,h=6)


### Genotyped SNPs
file.maf_distribution.genotyped <- "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/frq_files/Prote_370k_251011.genotyped_snps.no_mixup.with_ETypes.frq.cut5.txt"
df.maf_distribution.genotyped <- read.table(file.maf_distribution.genotyped, h=T)

p <- ggplot(df.maf_distribution.genotyped, aes(x=MAF)) + geom_histogram(binwidth=0.01)
p <- p + labs(title=sprintf("EGCUT genotyped SNPs. N_SNPs=%s", nrow(df.maf_distribution.genotyped)))
p
ggsave(file="maf_distribution_EGCUT_genotyped.pdf" ,w=8,h=6)
