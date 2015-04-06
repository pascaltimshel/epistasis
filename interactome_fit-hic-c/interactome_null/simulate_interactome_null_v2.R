############### SYNOPSIS ###################
# CLASS: script; data/table generative
# PURPOSE: Generate NULL table


###########################################

library(ggplot2)
library(permute)
library(tools) # for file_path_sans_ext()
library(dplyr)
library(reshape2)

rm(list=ls())

wd <- path.expand('~/git/epistasis/interactome_fit-hic-c/interactome_null')
setwd(wd)


###################### SOURCE ####################
source(file="function_generate_null_v2.R")
source(file="function_generate_null_utils_v2.R")
##################################################


###################### Read interaction table ####################
n_perm <- 1000
hic_cell_type <- "hESC" # "hESC" OR "hIMR90"
#hic_cell_type <- "hIMR90" # "hESC" OR "hIMR90"

### hIMR90
#p.q.threshold <- 1e-06 --> ??
#p.q.threshold <- 1e-07 --> 5 min
#p.q.threshold <- 1e-08 --> 1.2 min
#p.q.threshold <- 1e-09
#p.q.threshold <- 1e-10

### hESC
# 1e-13 --> XX min

list.timing.loop <- list()

for (p.q.threshold in c(1e-12,1e-13,1e-14,1e-15,1e-16,1e-17,1e-18)) { # hESC
#for (p.q.threshold in c(1e-13,1e-14,1e-15,1e-16,1e-17,1e-18)) { # hESC
#for (p.q.threshold in c(1e-16)) {
#for (p.q.threshold in c(1e-06, 1e-07, 1e-08, 1e-09, 1e-10)) { # IMR90
  time_start <- proc.time()
  
  str.path <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/tables_v1/interation_table.fit-hi-c.nosex.interchromosomal.%s.q_%s.txt" # e.g. interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-07.txt
  file.interaction_table <- sprintf(str.path, hic_cell_type, p.q.threshold)
  file.interaction_table
  
  df.interaction_table <- read.table(file=file.interaction_table, h=T, sep="\t", stringsAsFactors=F)
  
  
  ### Safety check: make sure the correct SORTED table is loaded ###
  stopifnot(all(df.interaction_table$chrA <= df.interaction_table$chrB)) # --> *MUST* be TRUE!
  
  
  ################################################################################################
  ###################################### Run function ########################################
  ################################################################################################
  r.list <- perm_restricted(n_perm, df.interaction_table)
  df.perm <- r.list[[1]]
  df.perm.stats <- r.list[[2]]
  ################################################################################################
  
  ################################################################################################
  ###################################### Make diagnostics ########################################
  ################################################################################################
  path.out.diagnostics = path.expand("~/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/interaction_tables/null_v2") # no trailing slash
  source(file="script_null_diagnostics_v2.R") # this script will run all the nessesary steps
  
  ################################################################################################
  ############################################ EXPORT ############################################
  ################################################################################################
  path.out = path.expand("~/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/interaction_tables/null_v2")
  write_df.perm(path.out=path.out)
  
  
  #####################################

  time_elapsed = proc.time() - time_start
  print(sprintf("Time elapsed, q=%s: %s", p.q.threshold, time_elapsed[3]))
  
  list.timing.loop[[as.character(p.q.threshold)]] <- time_elapsed[3]
}

list.timing.loop



