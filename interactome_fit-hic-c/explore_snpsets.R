############### SYNOPSIS ###################
# CLASS: script; explorative
# PURPOSE: read snpset table; estimate bonferroni correction
############################################

library(reshape2)
library(ggplot2)
library(plyr)

rm(list=ls())

wd <- "/Users/pascaltimshel/git/epistasis/interactome_fit-hic-c"
setwd(wd)
######################

################################# LOAD snpset table ############################


################################# LOADING interaction table ############################
#q.threshold <- 1e-6 # --> 26325
#q.threshold <- 1e-7 # --> 8114
#q.threshold <- 1e-8 # --> 3468
#q.threshold <- 1e-9 # --> 1021
#q.threshold <- 1e-10 # --> 432
#q.threshold <- 1e-11 # --> 214

str.path <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_snpsets/%s_snppool_hIMR90_q_%s/stats/df_interaction_table_snps.txt"
p.width <- 10000
p.q.threshold <- 1e-10
file.snpset_table <- sprintf(str.path, p.width, p.q.threshold)
file.snpset_table

df.snpset_table <- read.table(file.snpset_table, sep="\t", h=T, stringsAsFactors=F)
str(df.snpset_table)

######################### ESTIMATING bonferroni correcting #######################
n_tests <- df.snpset_table$setA_size * df.snpset_table$setB_size
n_tests.per_probe <- sum(n_tests) # this is the number of tests per probe
n_tests.per_probe
time.per_probe <- n_tests.per_probe/120000
time.per_probe
n_probes <- 10000 # number of probes
n_tests_total <- n_tests.per_probe * n_probes
n_tests_total

pval.nominal <- 0.05
pval.corrected <- pval.nominal/n_tests_total
pval.corrected

######################### P_values corrected #######################

### 1e-6 ###
# n_tests.per_probe = 6214029
# pval.corrected = 8.046309e-13

### 1e-10 ###
# n_tests.per_probe = 104492
# n_tests_total = 1044920000
# pval.corrected = 4.785055e-11


# 120 000 epistasis tests per second

######################### FUNCTION: threshold on q-value #######################
                                  # * DELETED *
#################################################################################


