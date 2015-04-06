############### SYNOPSIS ###################
# CLASS: script
# PURPOSE: sort and export fit-hi-c interaction table
############################################

library(reshape2)
library(ggplot2)
#library(plyr)
library(dplyr)

rm(list=ls())

wd <- "/Users/pascaltimshel/git/epistasis/interactome_fit-hic-c"
setwd(wd)
######################
################################# READ interaction table ############################

### Lan et al., inter-chromosomal, 
file.interaction_table <- "/Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Lan_et_al_chromosomal_interactions/lift_findItersection.intersection.paste.clean.nosex.updatedIDs.interchromosomal.fit-hi-c"
df.interaction_table <- read.table(file.interaction_table, sep="\t", h=T, stringsAsFactors=F)

### *CHANGE COLUMN NAMES*
colnames(df.interaction_table) <- c("chrA", "posA", "chrB", "posB", "contactCount", "p.value", "q.value")
head(df.interaction_table) # 9 columns
str(df.interaction_table)

### Check data integrity
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
  ## --> 0 removed for Lan et al. data

######################### QC data - duplicates removal ################

#x <- df.interaction_table %>% arrange(posA) # for manual inspection
#bool.duplicates <- duplicated(subset(df.interaction_table, select = -interactionID))
#sum(bool.duplicates)
  # --> NO duplicates

######################### QC data - outlier removal/"winsorizing" ################
### QC/outlier removal
# 1) we identify the "top hotspots": chromosomal coordinates with the most interactions (top X % most interactions)
# 2) we remove any interactions containing the a "top hotspot".
# NB: we use chr:pos as an identifier for the chrosmosomal coordinate.

### Adding TEMPORARY "index keys" columns to df.interaction_table
df.interaction_table$keyA <- with(df.interaction_table, paste0(chrA, ":", posA))
df.interaction_table$keyB <- with(df.interaction_table, paste0(chrB, ":", posB))

### Creating df with keyA and keyB VERTICALLY concattenated
df.vertcat <- data.frame(key = c(df.interaction_table$keyA, df.interaction_table$keyB))
nrow(df.vertcat) == 2*nrow(df.interaction_table) # --> TRUE

### df.key.table: counts for each key (interaction partner/chromosomal coordinate)
df.key.count <- as.data.frame(table(df.vertcat))
df.key.count <- df.key.count %>% arrange(desc(Freq))
df.key.count$percentage <- df.key.count$Freq/sum(df.key.count$Freq)*100

### get percentile/quantiles
probs <- seq(0,1,0.001) # steps of 0.1% (0.1/100)
df.quantile <- data.frame(quantile=probs, quantile_value=quantile(df.key.count$Freq, probs=probs))
df.quantile <- df.quantile %>% arrange(desc(quantile_value))
cutoff.count <- df.quantile[df.quantile$quantile==0.999, "quantile_value"]
cutoff.count

### get keys (chromosomal coordinates) that are MORE EXTREME (>) the cutoff (99.9%)
keys.outliers <- as.character(df.key.count[df.key.count$Freq > cutoff.count, 1])
length(keys.outliers) # --> 3
# df.vertcat	Freq
# 1	chr5:133628145	478
# 2	chr5:133628590	324
# 3	chr5:133625948	174
keys.outliers

### Removing outliers - subsetting dataframe
df.interaction_table.removed <- subset(df.interaction_table, keyA %in% keys.outliers | keyB %in% keys.outliers)
df.interaction_table.keep <- subset(df.interaction_table, !(keyA %in% keys.outliers | keyB %in% keys.outliers))
nrow(df.interaction_table.removed) # --> 976
nrow(df.interaction_table.keep) # --> 2521
nrow(df.interaction_table) == nrow(df.interaction_table.removed) + nrow(df.interaction_table.keep) # --> TRUE


#################################################################################

######################### MAIN LOOP - sorting table #######################
p.q.threshold <- c("OUTLIER_RM")
#q.threshold <- "OUTLIER_RM" # for manual "execution" of for loop
hic_cell_type <- "lan-et-al_K562"
df.interaction_table.input2mainLoop <- subset(df.interaction_table.keep, select = c(-keyA, -keyB))

for (q.threshold in p.q.threshold) {
  time_start <- proc.time()
  #df.interaction_table.sub.q <- df.interaction_table[df.interaction_table$q.value<=q.threshold,]
  df.interaction_table.sub.q <- df.interaction_table.input2mainLoop # ADJUSTED FOR Lan et al.
  n_int <- nrow(df.interaction_table.sub.q)
  cat(sprintf("q.threshold\t%s\t%s\n", q.threshold, n_int))
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
  stopifnot(with(df.loci.sort, all(chrA<=chrB)))
  
  ### Check that the positions are corectly swapped in the sorted data frame
  df.interaction_table.sub.q[316,]
  df.loci.sort[316,]
  
  ### *IMPORTANT*: replacing columns with sorted columns
  df.interaction_table.sub.q.sorted <- df.interaction_table.sub.q # copy
  df.interaction_table.sub.q.sorted[, c("chrA", "posA", "chrB", "posB")] <- df.loci.sort
  
  ######################### Check for duplicate interactions #############################
  ### Identifiying duplicates - first occurrence of record will *NOT* be marked as duplicate
  bool.duplicates <- duplicated(subset(df.interaction_table.sub.q.sorted, select = c("chrA", "posA", "chrB", "posB")))
  sum(bool.duplicates) # Lan et al --> 1258
  if (sum(bool.duplicates) != 0) {
    print(sprintf("Warning: detected %s duplicates in the interaction table", sum(bool.duplicates)))
    cat("Press [enter] to continue and duplicates will be removed")
    Sys.sleep(4)
    line <- readline()
  }
  
  ### Removing duplicates
  df.interaction_table.sub.q.sorted <- subset(df.interaction_table.sub.q.sorted, !bool.duplicates)
  print(sprintf("number of interactions AFTER duplicate removal: %s", nrow(df.interaction_table.sub.q.sorted)))
  
  ######################### EXPORTING interaction table #############################
  
  ### Adding UPDATED interactionID column
  df.interaction_table.sub.q.sorted$interactionID_final <- paste0("interaction_", 1:nrow(df.interaction_table.sub.q.sorted))
  str(df.interaction_table.sub.q.sorted)
  
  path.out <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/" # remember trailing backslash
  file <- sprintf("interation_table.fit-hi-c.nosex.interchromosomal.%s.q_%s.txt", hic_cell_type, q.threshold) # e.g. interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-06.txt
  file.out <- paste0(path.out, file)
  write.table(df.interaction_table.sub.q.sorted, file=file.out, sep="\t", quote=F, row.names=F, col.names=TRUE)

  time_elapsed = proc.time() - time_start
  print(sprintf("Time elapsed, q=%s: %s", q.threshold, time_elapsed[3]))
}

#str(df.loci.sort)
#str(df.interaction_table.sub.q.sorted)

##### Finding the counts/occurences for each genomic locations:
# a <- with(df.interaction_table, paste0(chrA, ":", posA))
# b <- with(df.interaction_table, paste0(chrB, ":", posB))
# x <- data.frame(t=c(a, b))
# x2 <- x %>% distinct() # --> 2021
# 
# df.table <- as.data.frame(table(x))
# df.table.s <- df.table %>% arrange(desc(Freq))


