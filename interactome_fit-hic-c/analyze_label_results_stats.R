############### SYNOPSIS ###################
# CLASS: script; summary stats
# PURPOSE: Analyze difference in "labelled" FastEpistasis epistatic pairs for Null and HiC.

# STEPS
# 1) Read in HiC labelled experiment/file
# 2) Read in "n_perm" (e.g. 1000) number of Null labelled experiments/files
# [POTENTIAL memory reduction solution: read in files in a loop and only keep the summary stats
# 3) Calculate summary statistics for each of these experiments/files
#   a) number of unique probes
#   b) {mean, variance, min, max} of epistatic SNP-pairs *per probe* 
#   c) {mean, variance, min, max} of epistatic SNP-pairs *per experiment interaction identifer*
#   d) [POTENTIALLY: plot distributions of these populations]
#   ALTERNATIVES: number of unique probes per *per experiment interaction identifer*
#   x) {mean/median, variance, min, max} of p-values
# 3) (make a few plots)

### THIGNS TO DO IN OTHER PROGRAMS ###
# Read in the *__results.epi.qt.lm.combined.txt file and look at the distribution of epistatic SNP-pairs over probes.

###########################################

library(ggplot2)


rm(list=ls())

wd <- '/cvar/jhlab/timshel/git/epistasis/interactome_fit-hic-c/interactome_null'
set(wd)
