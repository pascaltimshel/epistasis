############### SYNOPSIS ###################
# CLASS: THESIS PLOTS
# PURPOSE: Illustrate EMPIRICAL NULL
# Plots: density/histogram
############################################

library(ggplot2)
library(reshape2)

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

### Generate null distribution ()
p<-seq(0,1,0.01) # USED FOR THESIS: p<-seq(0,1,0.00001)
### F-distribution
df_numerator <- 5
density.f <- qf(p,df_numerator,1e3) + 100
density.norm <- qnorm(p, 100, 3)
df <- data.frame(p, density.f, density.norm)

p <- ggplot()
p <- p + geom_density(data=df, aes(x=density.f), alpha=1) # DENSITY | F-distribution
p <- p + geom_histogram(data=df, aes(x=density.f, y = ..density..), binwidth=0.1, fill="#1f78b4", alpha=0.8) # HISTOGRAM | F-distribution
#p <- p + geom_density(data=df, aes(x=density.norm, fill="norm"), alpha=0.5) # Normal
p <- p + geom_vline(xintercept = 102, linetype = "longdash", size=1, color="gray")
p <- p + geom_text(aes(x=102.5, y=0.5, label="s[obs]"), parse=TRUE) # OBS: parse=TRUE
p <- p + coord_cartesian(xlim=c(99.9, 104.5), ylim=c(0, 0.8)) 
p <- p + labs(x="Score", y="Frequency")
p
#ggsave(file="empirical_null_concept_histogram.pdf", w=8, h=4)

#####################################################
# Calculate p-value
p<-seq(0,1,0.01)
df_numerator <- 5
density.f <- qf(p,df_numerator,1e3)
sum(density.f >= 2)/length(density.f) # approximation: 0.07622924

pf(2, df_numerator,1e3, lower.tail=FALSE) # 0.0762263


#####################################################

### F-distibution
t <- qf(seq(0,1,0.00001), 6, 1e3)+120
qplot(t, geom="density", fill="blue")
### Normal distibution
t <- qnorm(seq(0,1,0.00001), 120,20)
qplot(t, geom="density", fill="blue")



#ggsave(file="empirical_null_concept_test.pdf", w=8, h=6)

######################################### stat_function METHOD - not so good ######################################
#ggplot(data.frame(x = c(-5, 5)), aes(x)) + stat_function(color="blue", fun = dnorm, args = list(mean = 2, sd = .5))

######################################### BARPLOT METHOD ######################################
### Generate null distribution ()
x<-seq(0,3,0.01)
df_numerator <- 5
hf.big <- df(x,df_numerator,10e5)
hchisq <- dchisq(x,df_numerator)
df.bar <- data.frame(x, hf.big, hchisq)

p <- ggplot()
p <- p + geom_bar(data=df.bar, aes(x=x, y=hf.big, fill="F(5,big)"), stat="identity", alpha=1)
p <- p + geom_density(data=df.bar, aes(x=x, y=hchisq, fill="hchisq"), stat="identity", alpha=1)

p


