############### SYNOPSIS ###################
# CLASS: script; quick; thesis stats
# PURPOSE: Count the number of unique interactions in each interaction table.
# FINAL GOAL: Calculate the genome coverage of each "parametrization".

#####################

############################################

library(reshape2)
library(ggplot2)
library(dplyr)
#library(plyr)
library(tools)

rm(list=ls())

wd <- path.expand("~/git/epistasis/thesis_stuff")
setwd(wd)
######################


path.input <- "/Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/interaction_tables/tables_v1"
files.interaction.tables <- list.files(path.input, full.names=TRUE)
files.interaction.tables

for (file in files.interaction.tables) {
  #file <- files.interaction.tables[1]
  file.basename <- basename(file)
  #print(file.basename)
  df <- read.table(file, h=T)
  interactions.A <- paste0(df$chrA,":",df$posA)
  interactions.B <- paste0(df$chrB,":",df$posB)
  interactions <- c(interactions.A, interactions.B)
  interactions.unique <- unique(interactions)
  cat(sprintf("%s\t%s\n", file.basename, length(interactions.unique)))
}

