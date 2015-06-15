############### SYNOPSIS ###################
# CLASS: script; analysis
# PURPOSE: run simulations for power detections

############################################

library(reshape2)
library(ggplot2)
library(dplyr)
library(tools)

rm(list=ls())

wd <- path.expand("~/git/epistasis/simulations")
setwd(wd)
############################################

