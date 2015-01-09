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
# --> 
file.fastEpi_table <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled/fastEpi_1e-10_0108_212539_cat_0109_121438/results.epi.qt.lm.combined"
df.fastEpi <- read.delim(file.fastEpi_table, stringsAsFactors=F)#
str(df.fastEpi)
nrow(df.fastEpi)

length(unique(df.fastEpi$SNP_A))

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

df.fastEpi.sorted_rsID.dup[1:2,]

### **** ###
df.fastEpi.sorted_rsID.unique <- unique(df.fastEpi.sorted_rsID)
df.fastEpi.sorted_rsID.unique

##############################################################################

##############################################################################

### Sort by p-value ###
df.fastEpi.sorted <- df.fastEpi[order(df.fastEpi$PVALUE), ]
head(df.fastEpi.sorted)

### Order of CHR: not sorted. This makes sense since we pooled all the SNPs from sorted interactions.
all(df.fastEpi$CHR <= df.fastEpi$CHR.1)
df.fastEpi$CHR <= df.fastEpi$CHR.1



#################################### GARBAGE ####################################

mat <- matrix(1:16, ncol=4)
mat
length(mat)

########### THIS WORKED - Extracting elements ##############
mat <- as.matrix(df.fastEpi.sorted_rsID[1:10, c("SNP_A", "SNP_B")])
mat
idx <- cbind(1:10, x.rsID.order[1:10,1])
idx
mat[idx]

SNP_A <- df.fastEpi.sorted_rsID[1:10, c("SNP_A", "SNP_B")][cbind(seq_along(x.rsID.order[1:10,1]), x.rsID.order[1:10,1])]
SNP_A

########### FAIL SORTING ATTEMPTS ##############
### fail attempt using apply()
apply(as.matrix(df.fastEpi.sorted_rsID[1:10, c("SNP_A", "SNP_B")]), 1, '[', x.rsID.order[1:10,1])

SNP_A <- df.fastEpi.sorted_rsID[1:10, c("SNP_A", "SNP_B")][[x.rsID.order[1:10,1]]]
SNP_A
df.fastEpi.sorted_rsID[1, c("SNP_A", "SNP_B")]
df.fastEpi.sorted_rsID[1, c("SNP_A", "SNP_B")][cbind(1,2)]

