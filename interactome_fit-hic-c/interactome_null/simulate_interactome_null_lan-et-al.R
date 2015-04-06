############### SYNOPSIS ###################
## COPY OF "simulate_interactome_null.R"
## This script if meant to process the interaction table from Lan et al

### NOTES
# See Evernote for a description of the criteria to the null.


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


###################### PARAMS ####################
p.q.threshold <- c("OUTLIER_RM")
hic_cell_type <- "lan-et-al_K562"
n_perm <- 1000

#################################### SCRIPT #################################

######### "FOR LOOP" - START #########
time_start <- proc.time()

str.path <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/tables_v1/interation_table.fit-hi-c.nosex.interchromosomal.%s.q_%s.txt" # e.g. interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-07.txt
file.interaction_table <- sprintf(str.path, hic_cell_type, p.q.threshold)
file.interaction_table

df.interaction_table <- read.table(file=file.interaction_table, h=T, sep="\t", stringsAsFactors=F)

### Safety check: make sure the correct SORTED table is loaded ###
stopifnot(all(df.interaction_table$chrA <= df.interaction_table$chrB)) # --> *MUST* be TRUE!

#################### Run function ####################
r.list <- perm_restricted(n_perm, df.interaction_table)
df.perm <- r.list[[1]]
df.perm.stats <- r.list[[2]]
######################################################


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


time_elapsed = proc.time() - time_start
print(sprintf("Time elapsed, q=%s: %s", p.q.threshold, time_elapsed[3]))

######### "FOR LOOP" - END #########
