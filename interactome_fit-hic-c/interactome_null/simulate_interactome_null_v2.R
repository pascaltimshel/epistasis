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

################################################################################################
###################################### Broad/OSX switch ########################################
################################################################################################

#os.execution <- "broad"
os.execution <- "osx"

# Broad NOTES
# reuse R-3.1
# bsub -J simulate_interactome_null_v2-broad-hESC -o simulate_interactome_null_v2-broad-hESC.bsub.out -q MEDPOP -n 1 -R 'rusage[mem=10]' Rscript simulate_interactome_null_v2.R

###################################### Switch ###################################### 
if (os.execution == "osx") {
  wd <- path.expand('~/git/epistasis/interactome_fit-hic-c/interactome_null')
  setwd(wd)
  
  ### Set paths
  path.interaction_tables <- path.expand("~/p_HiC/Ferhat_Ay_2014/interaction_tables/tables_v1")
  path.out.diagnostics.base = path.expand("~/p_HiC/Ferhat_Ay_2014/interaction_tables/null_v2")
  path.out.export = path.expand("~/p_HiC/Ferhat_Ay_2014/interaction_tables/null_v2")
  #path.out.diagnostics.base = path.expand("~/p_HiC/Ferhat_Ay_2014/interaction_tables/null_v2_small")
  #path.out.export = path.expand("~/p_HiC/Ferhat_Ay_2014/interaction_tables/null_v2_small")
  
} else if (os.execution == "broad") {
  options(echo=TRUE)
  wd <- path.expand('/cvar/jhlab/timshel/git/epistasis/interactome_fit-hic-c/interactome_null')
  setwd(wd)
  
  ### Set paths
  path.interaction_tables <- path.expand("/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/interaction_tables/tables_v1")
  path.out.diagnostics.base = path.expand("/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/interaction_tables/null_v2")
  path.out.export = path.expand("/cvar/jhlab/timshel/egcut/interactome_fit-hi-c/interaction_tables/null_v2")
  
} else {
  stop("Did not recognize os.execution argument")
}
################################################################################################
################################################################################################


###################### SOURCE ####################
#source(file="function_generate_null_v2_parallel.R")
source(file="function_generate_null_v2.R")
source(file="function_generate_null_utils_v2.R")
##################################################


###################### Set params ####################
n_perm <- 1000
#n_perm <- 100

######### Set cell type | this is only used for the OUTPUT filename #########
#hic_cell_type <- "hESC" # q_lt_0.001.inter | "hESC" "hIMR90"
#hic_cell_type <- "hESC-contactCount_1" # contactCount_1
hic_cell_type <- "hIMR90-contactCount_1" # contactCount_1

### hIMR90
# 1e-06 [26 k interactions] --> 70 min
# 1e-07 --> 5 min
# 1e-08 --> 1.2 min
# 1e-09
# 1e-10

### hESC
# 1e-12 [39 k interactions] --> 192 min (3.2 h)
# 1e-13 [~25 k interactions] --> 56 min
# 1e-14 --> 10 min

### contactCount_1
# hESC head_5001 [5k interaction] --> 13 min
# hIMR90 head_5001 [5k interaction] --> 6 min

list.timing.loop <- list()

#for (p.q.threshold in c(1e-12,1e-13,1e-14,1e-15,1e-16,1e-17,1e-18)) { # hESC
#for (p.q.threshold in c(1e-13,1e-14,1e-15,1e-16,1e-17,1e-18)) { # hESC
#for (p.q.threshold in c(1e-17,1e-18)) {
#for (p.q.threshold in c(1e-06, 1e-07, 1e-08, 1e-09, 1e-10)) { # IMR90
#for (p.q.threshold in c(1e-08, 1e-09, 1e-10)) { # IMR90
for (p.q.threshold in c(1)) { # contactCount_1 | make sure that the given "q.threshold" exists
  time_start <- proc.time()
  
  #str.path <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/tables_v1/interation_table.fit-hi-c.nosex.interchromosomal.%s.q_%s.txt" # e.g. interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-07.txt
  str.path <- file.path(path.interaction_tables, "interation_table.fit-hi-c.nosex.interchromosomal.%s.q_%s.txt") # e.g. interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-07.txt
  file.interaction_table <- sprintf(str.path, hic_cell_type, p.q.threshold)
  print(file.interaction_table)
  
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
  # REMEMBER TO SET "path.out.diagnostics" BEFORE SOURCING/RUNNING "script_null_diagnostics_v2.R"
  path.out.diagnostics <- file.path(path.out.diagnostics.base, hic_cell_type)
  if (file.exists(path.out.diagnostics)) {
    print("OBS: path.out.diagnostics folder exists")
  }
  dir.create(path.out.diagnostics, showWarnings=FALSE) # dir.create() does NOT crash if the directory already exists, it just prints out a warning.
  print(path.out.diagnostics) # no trailing slash
  
  # this script will run all the nessesary steps for generating DIAGNOSTIC PLOTS
  source(file="script_null_diagnostics_v2.R")
  
  ################################################################################################
  ############################################ EXPORT ############################################
  ################################################################################################
  ### REMEMBER to set path.out.export before running function
  write_df.perm(path.out=path.out.export)
  
  
  #####################################

  time_elapsed = proc.time() - time_start
  print(sprintf("Time elapsed, q=%s: %s", p.q.threshold, time_elapsed[3]))
  
  list.timing.loop[[as.character(p.q.threshold)]] <- time_elapsed[3]
}

print(list.timing.loop)



