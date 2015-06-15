############### SYNOPSIS ###################
# CLASS: THESIS PLOTS
# PURPOSE: Plot MANY histogram plots for empirical p-values of spatial epistasis
# Plots: Histograms
############################################

library(ggplot2)
library(reshape2)
library(plyr) # for rename()

#########
rm(list=ls())

wd <- "~/git/epistasis/thesis_plots"
setwd(wd)

#########
################## XXX ###################

### Set file name
run.epi <- "hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8"
file.in <- sprintf("/Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/compiled_2015-06-15_epistatic_counts_csv/epistatic_counts_%s.csv", run.epi)
### Read file
df <- read.csv(file.in)
str(df)
### Rename variables
rename(df, c("X"="EID", 
             "count_significant"="No filter",
             "count_significant_pruned_EIID"="Filter 1 (EIID)",
             "count_significant_pruned_hemani"="Filter 2 (Hemani)"
             )
       )
### Melt
df.melt <- melt(df) # value, variable
### REORDER VARIABLES to define the order of the facet_wrap()
df.melt$variable <- factor(df.melt$variable, levels = c("small", "medium", "large")) #  could also use relevel()
str(df.melt)
### Plot
p <- ggplot(df.melt)
p <- p + geom_histogram(aes(x=value), binwidth=1, alpha=1)
p <- p + facet_wrap(~variable, scales="free", ncol=1)
# To display different lines in different facets, you need to create a data frame.
p <- p + geom_vline(aes(xintercept=)

p

plot.w <- 8
plot.h <- 4
plot.filename <- sprintf("epistasis_enrichment_%s-%s-%s.pdf", run.epi, plot.w, plot.h)
ggsave(file=plot.filename, w=plot.w, h=plot.h)


