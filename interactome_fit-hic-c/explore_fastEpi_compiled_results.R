############### SYNOPSIS ###################
# CLASS: script; explorative
# PURPOSE: read fastEpi compiled table
# x) number of unique probes and snps. 
# x) any duplicated probe/SNP pair?

#snpStats: Replicate p-value using F-test
#snpStats: Does SNPs that are often detected as significant have a high missingness rate?
#snpStats: Are the SNPs in LD?
############################################

library(reshape2)
library(ggplot2)
library(plyr)

rm(list=ls())

wd <- "/Users/pascaltimshel/git/epistasis/interactome_fit-hic-c"
setwd(wd)
######################

################################# LOAD fastEpi table ############################
### EXAMPLE ###
#CHR  SNP_A  CHR	SNP_B	BETA	CHISQ	PVALUE	PHENOTYPE
#9	rs4579584	20	rs2426193	-0.25225	49.39828	2.08928E-12	ILMN_1799744
#9	rs4579584	20	rs7266126	-0.25303	49.61180	1.87383E-12	ILMN_1799744
#4	rs7681358	11	rs17137547	-1.49101	72.67622	1.52763E-17	ILMN_1694385

##################### Initial TABLES #####################
### COMPLETE TABLE
#file.fastEpi_table <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/fastEpi_1e-10_0108_212539/results.epi.qt.lm.combined.stripped"
  # --> nrow=833144

##################### Current TABLES #####################
file.fastEpi_table <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8/hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8__150122_113919___results.epi.qt.lm.combined.txt"
  # --> nrow=49827

### READ table
df.fastEpi <- read.delim(file.fastEpi_table, stringsAsFactors=F)#
str(df.fastEpi)
nrow(df.fastEpi)

### UNIQUE SNPs and probes
SNP_A.unique <- unique(df.fastEpi$SNP_A); length(SNP_A.unique) # --> 6512 | new=2083
SNP_B.unique <- unique(df.fastEpi$SNP_B); length(SNP_B.unique) # --> 6501 | new=2026
SNP_AB.intersect <- intersect(SNP_A.unique, SNP_B.unique); length(SNP_AB.intersect) # --> 54 | new=10
SNP_AB.union <- union(SNP_A.unique, SNP_B.unique); length(SNP_AB.union) # --> 12959 | new=4090
probes.unique <- unique(df.fastEpi$PHENOTYPE); length(probes.unique) # --> 1127 | new=1224
# extract UNION snps
#cat(SNP_AB.union, sep="\n", file=paste0(dirname(file.fastEpi_table),"/SNP_AB.union.txt"))


### Count number of pairs below a certain p-value threshold
n.passing <- sum(df.fastEpi$PVALUE<1e-11)
n.passing/nrow(df.fastEpi)*100 # Percentage

### Checking for duplicates ###
# --> *We need a sorted rsID to do this*
anyDuplicated(df.fastEpi[,c("SNP_A", "SNP_B", "PHENOTYPE")]) # --> 0 duplicated. This is the "unsorted" data frame, so we cannot be sure.
x.rsID <- cbind(df.fastEpi$SNP_A, df.fastEpi$SNP_B) #x.rsID <- matrix(c(df.fastEpi$SNP_A, df.fastEpi$SNP_B), ncol=2) # gives the same as cbind()
head(x.rsID)
x.rsID.order <- t(apply(x.rsID, 1, order))
head(x.rsID.order)
### Copy data frame
df.fastEpi.sorted_rsID <- df.fastEpi # this data frame MUST have the same ROW order (keep them in sync)
### Creating extraction indices
# NOTE: the length of a matrix is the NUMBER OF ELEMENTS
idxA <- cbind(1:nrow(x.rsID.order), x.rsID.order[,1]) # first column of x.rsID.order
idxB <- cbind(1:nrow(x.rsID.order), x.rsID.order[,2]) # second column of x.rsID.order
### "Swapping" rsIDs: *OBS*: we are using df.fastEpi because after modifying "SNP_A" in df.fastEpi.sorted_rsID the indexing is no longer valid.
# it is important that the columns are extracted in the CORRECT ORDER c("SNP_A", "SNP_B")
df.fastEpi.sorted_rsID$SNP_A <- df.fastEpi[, c("SNP_A", "SNP_B")][idxA]
df.fastEpi.sorted_rsID$SNP_B <- df.fastEpi[, c("SNP_A", "SNP_B")][idxB]
### "Swapping" CHR: we are using df.fastEpi because after modifying "SNP_A" in df.fastEpi.sorted_rsID the indexing is no longer valid.
# it is important that the columns are extracted in the CORRECT ORDER c("CHR", "CHR.1")
df.fastEpi.sorted_rsID$CHR <- df.fastEpi[, c("CHR", "CHR.1")][idxA]
df.fastEpi.sorted_rsID$CHR.1 <- df.fastEpi[, c("CHR", "CHR.1")][idxB]
### Inspecting results
head(df.fastEpi.sorted_rsID)
head(df.fastEpi)
### Checking for duplicates again ###
# *OBS*: here we are only checking for DUPLICATED SNP-probe pairs. That is, the p-values does NOT have to be the same.
any(duplicated(df.fastEpi.sorted_rsID[,c("SNP_A", "SNP_B", "PHENOTYPE")])) # anyDuplicated() --> 14524
sum(duplicated(df.fastEpi.sorted_rsID[,c("SNP_A", "SNP_B", "PHENOTYPE")])) # --> 30
sum(duplicated(df.fastEpi.sorted_rsID)) #  all duplicated 
idx.duplicated <- duplicated(df.fastEpi.sorted_rsID[,c("SNP_A", "SNP_B", "PHENOTYPE")]) | duplicated(df.fastEpi.sorted_rsID[,c("SNP_A", "SNP_B", "PHENOTYPE")], fromLast=T)
df.fastEpi.sorted_rsID.dup <- df.fastEpi.sorted_rsID[idx.duplicated, ]
df.fastEpi.sorted_rsID.dup <- df.fastEpi.sorted_rsID.dup[with(df.fastEpi.sorted_rsID.dup, order(SNP_A, SNP_B)),]

### write csv
file.duplicated <- paste0(dirname(file.fastEpi_table), "/", "duplicated_SNP-probe_pairs.csv")
file.duplicated
#write.csv(df.fastEpi.sorted_rsID.dup, file=file.duplicated)


### **** ###
df.fastEpi.sorted_rsID.unique <- unique(df.fastEpi.sorted_rsID)
df.fastEpi.sorted_rsID.unique

##############################################################################

##############################################################################

### Sort by p-value ###
df.fastEpi.sorted <- df.fastEpi[order(df.fastEpi$PVALUE), ]
head(df.fastEpi.sorted)
### Write to csv
file.fastEpi_table.p_value.sorted <- paste0(dirname(file.fastEpi_table), "/", "fastEpi_table.p_value.sorted.csv")
file.fastEpi_table.p_value.sorted
n.top <- 10000
write.csv(df.fastEpi.sorted, file=file.fastEpi_table.p_value.sorted, row.names=F)
write.csv(df.fastEpi.sorted[1:n.top,], file=paste0(dirname(file.fastEpi_table), "/", "fastEpi_table.p_value.sorted.top1000.csv"), row.names=F)


n.top <- 1000
df.fastEpi.sorted.n_top <- df.fastEpi.sorted[1:n.top, ]
SNP_A.unique <- unique(df.fastEpi.sorted.n_top$SNP_A); length(SNP_A.unique) # --> 47
SNP_B.unique <- unique(df.fastEpi.sorted.n_top$SNP_B); length(SNP_B.unique) # --> 37
SNP_AB.intersect <- intersect(SNP_A.unique, SNP_B.unique); length(SNP_AB.intersect) # --> 0










#################################### OBSERVATIONS ####################################
### Order of CHR: not sorted. This makes sense since we pooled all the SNPs from sorted interactions.
all(df.fastEpi$CHR <= df.fastEpi$CHR.1)
df.fastEpi$CHR <= df.fastEpi$CHR.1
















# #################################### GARBAGE ####################################
# 
# mat <- matrix(1:16, ncol=4)
# mat
# length(mat)
# 
# ########### THIS WORKED - Extracting elements ##############
# mat <- as.matrix(df.fastEpi.sorted_rsID[1:10, c("SNP_A", "SNP_B")])
# mat
# idx <- cbind(1:10, x.rsID.order[1:10,1])
# idx
# mat[idx]
# 
# SNP_A <- df.fastEpi.sorted_rsID[1:10, c("SNP_A", "SNP_B")][cbind(seq_along(x.rsID.order[1:10,1]), x.rsID.order[1:10,1])]
# SNP_A
# 
# ########### FAIL SORTING ATTEMPTS ##############
# ### fail attempt using apply()
# apply(as.matrix(df.fastEpi.sorted_rsID[1:10, c("SNP_A", "SNP_B")]), 1, '[', x.rsID.order[1:10,1])
# 
# SNP_A <- df.fastEpi.sorted_rsID[1:10, c("SNP_A", "SNP_B")][[x.rsID.order[1:10,1]]]
# SNP_A
# df.fastEpi.sorted_rsID[1, c("SNP_A", "SNP_B")]
# df.fastEpi.sorted_rsID[1, c("SNP_A", "SNP_B")][cbind(1,2)]
# 
