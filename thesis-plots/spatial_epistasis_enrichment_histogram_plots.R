############### SYNOPSIS ###################
# CLASS: THESIS PLOTS
# PURPOSE: Plot MANY histogram plots for empirical p-values of spatial epistasis
# Plots: Histograms

### OBS: this script will have plots in the WD from setwd()
############################################

library(ggplot2)
library(reshape2)
library(plyr) # for rename()
library(dplyr)

#########
rm(list=ls())

wd <- "~/git/epistasis/thesis_plots"
setwd(wd)

#########
################## FUNCTION ###################
plot_histograms <- function(run.epi){
  print(sprintf("running: %s", run.epi))
  ### Set file name
  #run.epi <- "hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8"
  ### OBS: remember to UPDATE THE file.in path!
  file.in <- sprintf("/Users/pascaltimshel/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/compiled_2015-06-29_epistatic_counts_csv/epistatic_counts_%s.csv", run.epi)
  ### Read file
  df <- read.csv(file.in)
  #str(df)
  
  ## SUBSET df (if needed)
  #cols2keep <- c("X", "count_significant", "count_significant_pruned_EIID_probe", "count_significant_pruned_hemani")
  cols2keep <- c("X", "count_significant", "count_significant_pruned_EIID_probe")
  df <- subset(df, select=cols2keep)
  
  
  ### Rename variables
  df <- plyr::rename(df, c("X"="EID", 
               "count_significant"="No filter",
               #"count_significant_pruned_EIID"="Filter 1",
               "count_significant_pruned_EIID_probe"="Filter (interaction pair)"
               #"count_significant_pruned_hemani"="Hemani et al. filter (chromosome pair)"
               )
         )
  
  
  
  ### Melt
  df.melt <- melt(df, id.vars="EID") # value, variable
  
  ### REORDER VARIABLES to define the order of the facet_wrap()
  #df.melt$variable <- factor(df.melt$variable, levels = c("No filter", "Filter 1", "Filter 2")) #  could also use relevel()
  #df.melt$variable <- factor(df.melt$variable, levels = c("Without filter", "With filter")) #  could also use relevel()
  df.melt$variable <- factor(df.melt$variable, levels = c("No filter", "Filter (interaction pair)")) #  could also use relevel()
  #df.melt$variable <- factor(df.melt$variable, levels = c("No filter", "Filter (interaction pair)", "Hemani et al. filter (chromosome pair)")) #  could also use relevel()
  #str(df.melt)
  
  ### Extract hic_1 observations
  df.melt.hic <- subset(df.melt, EID=="hic_1")
  ### Exclude hic_1 observation (we do not want to show this observation in the histogram)
  df.melt <- subset(df.melt, EID!="hic_1")
  nrow(df.melt) == 1000*(ncol(df)-1) # MUST be TRUE
  
  ### Calculate null mean SNP-probe pairs (with hic_1 excluded)
  df.melt.null.mean <- df.melt %>% group_by(variable) %>% summarise(mean=mean(value), max_y=max(value))
  # Construct label to be used in ggplot geom_text with parse=TRUE
  #df.melt.null.mean$label <- sapply(df.melt.null.mean$mean, function(val){sprintf("mu==%s", round(val))}) # THIS WORKS
  df.melt.null.mean$label <- paste0("mu==",round(df.melt.null.mean$mean)) # THIS ALSO WORKS
  df.melt.null.mean
  
  ### Calculate max count for the null - TO BE USED FOR geom_text y-position
  
  
  ### Plot
  p <- ggplot(df.melt)
  p <- p + geom_histogram(aes(x=value), binwidth=1, alpha=1)
  p <- p + facet_wrap(~variable, scales="free", ncol=1)
  #### To display different lines in different facets, you need to create a data frame.
  ### Observed Hi-C count | vline + text
  p <- p + geom_vline(data=df.melt.hic, aes(xintercept=value), color="red")
  p <- p + geom_text(data=df.melt.hic, aes(x=value, y=0, label=value), hjust=1, vjust=0, color="red")
  ### Mean null count | vline + text
  # ggplot_build(p) *SEMI HACK*
  
  df.tmp.hack <- ggplot_build(p)$data[[1]][c("count", "PANEL", "group")] %>% group_by(PANEL) %>% summarise(max_count = max(count))
  df.hack <- df.tmp.hack %>% inner_join(ggplot_build(p)$panel$layout, by="PANEL")
  df.melt.null.mean.plot <- df.melt.null.mean %>% inner_join(df.hack, by="variable")
  p <- p + geom_vline(data=df.melt.null.mean.plot, aes(xintercept=mean), color="black", linetype="dashed", size=1.2)
  p <- p + geom_text(data=df.melt.null.mean.plot, aes(x=mean, y=max_count, label=label), parse=TRUE, hjust=1, vjust=1, color="white")
  
  ### Labels
  p <- p + labs(x="Number of epistatic SNP-probe pairs", y="Count")
  p
  
  
  
  plot.w <- 6
  plot.h <- 3
  plot.filename <- sprintf("epistasis_enrichment_%s-%sx%s.pdf", run.epi, plot.w, plot.h)
  ggsave(file=plot.filename, w=plot.w, h=plot.h)
}

################## RUN FUNCTION ###################
epi.runs <-  c("hESC_width_500_maf_5_q_1e-14_epi1_1e-8",
               "hESC_width_500_maf_5_q_1e-16_epi1_1e-8",
               "hESC_width_1000_maf_5_q_1e-12_epi1_1e-10",
               "hESC_width_2500_maf_5_q_1e-13_epi1_1e-10",
               "hESC-contactCount_1_width_1000_maf_5_q_1_epi1_1e-8",
               "hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8",
               "hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8",
               "hIMR90_width_1000_maf_5_q_1e-06_epi1_1e-10",
               "hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8",
               "hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10",
               "hIMR90-contactCount_1_width_1000_maf_5_q_1_epi1_1e-8",
               "lan-et-al_K562_width_1000_maf_5_q_OUTLIER_RM_epi1_1e-8",
               "lan-et-al_K562_width_5000_maf_5_q_OUTLIER_RM_epi1_1e-8")

stopifnot(length(epi.runs)==13)

### Run function
sapply(epi.runs, plot_histograms)

