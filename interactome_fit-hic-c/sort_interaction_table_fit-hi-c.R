############### SYNOPSIS ###################
# CLASS: script
# PURPOSE: sort and export fit-hi-c interaction table
############################################

library(reshape2)
library(ggplot2)
library(plyr)

rm(list=ls())

wd <- "/Users/pascaltimshel/git/epistasis/interactome_fit-hic-c"
setwd(wd)
######################

################################# STATS interaction table ############################
### hESC_HindIII - ~150 MB file

### hIMR90_HindIII - ~60 MB file
################################# LOADING interaction table ############################

### hESC_HindIII
#load(file="/Users/pascaltimshel/Dropbox/5_Data/p_HiC_viz/RData/hESC_HindIII_hg19.spline_pass1.res10000.significances.q_lt_0.001.inter.RData") # df.interaction_table.hESC

### hIMR90_HindIII
load(file="/Users/pascaltimshel/Dropbox/5_Data/p_HiC_viz/RData/hIMR90_HindIII_hg19.spline_pass1.res10000.significances.q_lt_0.001.inter.RData") # df.interaction_table.hIMR90

### Switch - THIS MUST MATCH WITH the loaded data
#df.interaction_table <- df.interaction_table.hESC
df.interaction_table <- df.interaction_table.hIMR90

### Check data integrity
anyNA(df.interaction_table) # False!
head(df.interaction_table) # 9 columns

### *CHANGE COLUMN NAMES*
colnames(df.interaction_table) <- c("chrA", "posA", "chrB", "posB", "contactCount", "p.value", "q.value")
head(df.interaction_table) # 9 columns
str(df.interaction_table)

### Adding interactionID column
df.interaction_table$interactionID <- paste0("interaction_", 1:nrow(df.interaction_table))
str(df.interaction_table)

######################### Removing chrX and chrY #######################

### NEW | *remove chrX and chrY*
bool.chrA.sex <- df.interaction_table$chrA %in% c("chrX", "chrY")
bool.chrB.sex <- df.interaction_table$chrB %in% c("chrX", "chrY")
sum(bool.chrA.sex); sum(bool.chrB.sex)
sum(bool.chrA.sex | bool.chrB.sex) # OR-logic
df.interaction_table <- df.interaction_table[!(bool.chrA.sex | bool.chrB.sex), ]


######################### FUNCTION: threshold on q-value #######################
                                  # * DELETED *
#################################################################################

######################### SUBSETTING: threshold on q-value #######################
q.threshold <- 1e-8 # --> 3468
#q.threshold <- 1e-9 # --> 1021
#q.threshold <- 1e-10 # --> 432
#q.threshold <- 1e-11 # --> 432
df.interaction_table.sub.q <- df.interaction_table[df.interaction_table$q.value<=q.threshold,]
nrow(df.interaction_table.sub.q)
head(df.interaction_table.sub.q)


######################### STRIPPING leading "chr" in chrA and chrB #############################
df.interaction_table.sub.q$chrA <- as.numeric(substring(df.interaction_table.sub.q$chrA, 4))
df.interaction_table.sub.q$chrB <- as.numeric(substring(df.interaction_table.sub.q$chrB, 4))
str(df.interaction_table.sub.q)


######################### Sorting interaction table #############################
source("function_sort_interaction_table_fit-hi-c.R")

### USING sort_interactions() - sort chrA and chrB bases on their numeric value
str(df.interaction_table.sub.q)
df.loci.sort <- sort_interactions(df.interaction_table.sub.q[, c("chrA", "posA", "chrB", "posB")])
### Check: chrA<=chrB
with(df.loci.sort, all(chrA<=chrB))

### Check that the positions are corectly swapped in the sorted data frame
df.interaction_table.sub.q[316,]
df.loci.sort[316,]

### *IMPORTANT*: replacing columns with sorted columns
df.interaction_table.sub.q.sorted <- df.interaction_table.sub.q # copy
df.interaction_table.sub.q.sorted[, c("chrA", "posA", "chrB", "posB")] <- df.loci.sort

######################### EXPORTING interaction table #############################

### Adding UPDATED interactionID column
df.interaction_table.sub.q.sorted$interactionID_final <- paste0("interaction_", 1:nrow(df.interaction_table.sub.q.sorted))
str(df.interaction_table.sub.q.sorted)

path.out <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/" # remember trailing backslash
file <- sprintf("interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_%s.txt", q.threshold)
file.out <- paste0(path.out, file)
write.table(df.interaction_table.sub.q.sorted, file=file.out, sep="\t", quote=F, row.names=F, col.names=TRUE)

