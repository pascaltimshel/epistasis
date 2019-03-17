############### SYNOPSIS ###################
# CLASS: THESIS PLOTS
# PURPOSE: Illustrate epistasis
# Plots: Line plot
############################################

library(ggplot2)
library(reshape2)

library(combinat)
#########
rm(list=ls())

wd <- "~/git/epistasis/thesis_plots"
setwd(wd)

#########
################## Generate data frame ###################

setup_plot <- function(p) {
  p <- p + labs(x="SNP 1", y="Phenotype")
  p <- p + guides(color=guide_legend(title="SNP 2"))
  #p <- p + theme_bw() # remove gray background
  return(p)
}

### Generate Allele Codes
allele.codes <- rep(c(0,1,2), 3)
df.alleles <- unique(expand.grid(snp1=allele.codes, snp2=allele.codes))

### Set beta coefficients - Only main effect from snp2 | NO EPISTASIS
beta.1 <- 0
beta.2 <- 0.3
beta.12 <- 0
### Set response variable (phenotype)
y1 <- with(df.alleles, beta.1*snp1 + beta.2*snp2 + beta.12*snp1*snp2)
df.epi <- with(df.alleles, data.frame(snp1=as.factor(snp1), snp2=as.factor(snp2), y=y1))
### Plot
p <- ggplot(df.epi, aes(x=snp1, group=snp2, color=snp2, y=y)) + geom_point() + geom_line()
p <- setup_plot(p)
ggsave(file="epistasis_definition_line_plot_main_snp2.pdf", w=6, h=3)



### Set beta coefficients - Only main effect from snp1 and snp2 | NO EPISTASIS
beta.1 <- 0.1
beta.2 <- 0.3
beta.12 <- 0
### Set response variable (phenotype)
y1 <- with(df.alleles, beta.1*snp1 + beta.2*snp2 + beta.12*snp1*snp2)
df.epi <- with(df.alleles, data.frame(snp1=as.factor(snp1), snp2=as.factor(snp2), y=y1))
### Plot
p <- ggplot(df.epi, aes(x=snp1, group=snp2, color=snp2, y=y)) + geom_point() + geom_line()
p <- setup_plot(p)
ggsave(file="epistasis_definition_line_plot_main_snp1_snp2.pdf", w=6, h=3)


### Set beta coefficients - Only main effect from snp1 and snp2 | NO EPISTASIS
beta.1 <- 0.1
beta.2 <- 0.3
beta.12 <- 0.5
### Set response variable (phenotype)
y1 <- with(df.alleles, beta.1*snp1 + beta.2*snp2 + beta.12*snp1*snp2)
df.epi <- with(df.alleles, data.frame(snp1=as.factor(snp1), snp2=as.factor(snp2), y=y1))
### Plot
p <- ggplot(df.epi, aes(x=snp1, group=snp2, color=snp2, y=y)) + geom_point() + geom_line()
p <- setup_plot(p)
p
ggsave(file="epistasis_definition_line_plot_main_snp1_snp2_snp12.pdf", w=6, h=3)


