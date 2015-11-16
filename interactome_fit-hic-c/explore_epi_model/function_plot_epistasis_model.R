############### SYNOPSIS ###################
# *IMPORTANT*: this script is a MODIFIED DUPLICATE of the "..XXX/EGCUT_DATA/R_epistasis_tests" directory

### Modifications
# 1) added argument: a data frame to extract snp1, snp2 and probename is parsed instead of using "df.snp_probe_pairs.subset"
# 2) *NOT YET IMPLEMENTED* added argument: directory to save plots
# 3) using "ggsave()" instead of "print(p)" dev.
# 5) added argument: "add_data_heatmap". Plot the data table with the number of missing values.
  # 5.1) sourcing ggplot "multiplot" function

# CLASS: function
# PURPOSE: plot epistasis model
# INPUT:
  # 1) df: a data frame with the columns: snp1, snp2, probename 
  # 2) idx: the range/indices of observations to plot from the data frame
  # 3) save_images and save_significant_treshold: determines whether to save the plot or not.
# OUTPUT: optionally: saved pdf plots in the *CURRENT DIRECTORY*.

# Layers in plot
# 1) points: expression values for individuals
# 2) line: mean of expression values
# 3) dashed line: fitted epistasis model
# Finaly, a grid is wrapped around it
############################################

plot_EpiModel <- function(idx, df.input, path.out="./", add_text_annotation=TRUE, add_data_heatmap=FALSE, save_images=TRUE, save_significant_treshold=1, plot.points.alpha=0.3) {
  ### Defaults
  # path.out --> current working directory.
  
  
  #plots <- list()
  
  ### FOR DEBUGGING | "STABLE" CASE ##
  #i <- 2
  #df.input <- df.fastEpistasis.results
  
  ### FOR DEBUGGING | TEMPORARY CASE ##
  i <- 2
  df.input <- df.selected
  
  for (i in idx) {
    print(sprintf("plotting #%d/#%d", i, length(idx)))
    
    name.probe <- df.input[i, "probename"]
    name.snp1 <-df.input[i, "snp1"]
    name.snp2 <- df.input[i, "snp2"]
    
    
    ### OLD APPROACH: no parsing of a data frame argument
    #name.probe <- df.snp_probe_pairs.subset[i, "probename"]
    #name.snp1 <-df.snp_probe_pairs.subset[i, "snp1"]
    #name.snp2 <- df.snp_probe_pairs.subset[i, "snp2"]
    
    # Extract data
    snp1 <- geno[, colnames(geno) == name.snp1]
    snp2 <- geno[, colnames(geno) == name.snp2]
    probe <- df.probes[, colnames(df.probes) == name.probe]
    
    print(name.snp1)
    print(name.snp2)
    
    ### Misc
    #tab <- table(data.frame(snp1, snp2))
    #tab
    snp_ld <- cor(snp1, snp2, use="pairwise.complete.obs")
    n_usable_data_points <- sum(!is.na(snp1) & !is.na(snp2))
    n_missing_data_points <- sum(is.na(snp1)) + sum(is.na(snp2))
    
    print(i)
    
    #### Multiplicative model
    fullmod.multiplicative <- lm(probe ~ snp1 + snp2 + snp1:snp2) # Full
    margmod.multiplicative <- lm(probe ~ snp1 + snp2) # Reduced
    fullmod.multiplicative
    summary(fullmod.multiplicative)
    summary(fullmod.multiplicative)$r.squared # Multiple R-squared
    beta_mu <- fullmod.multiplicative$coefficients[1]
    beta_snp1 <- fullmod.multiplicative$coefficients[2]
    beta_snp2 <- fullmod.multiplicative$coefficients[3]
    beta_int <- fullmod.multiplicative$coefficients[4] # coefficient for interaction term (snp1:snp2)
    inttest_multiplicative <- anova(margmod.multiplicative, fullmod.multiplicative) # Test for model reduction
    inttest_multiplicative_pval <- inttest_multiplicative$P[2]
    inttest_multiplicative_pval
    #-log10(inttest_multiplicative_pval) # P-value
    
    #leveragePlots(fullmod.multiplicative)
    #plot(fullmod.multiplicative)
    
    ### Summarizing expression data
    df <- data.frame(snp1=as.factor(snp1), snp2=as.factor(snp2), probe=probe)
    df.mean <- ddply(df, .(snp1, snp2), summarize,
                     mean=mean(probe, na.rm=T),
                     sd=sd(probe, na.rm=T))
    ### Replace NA in sd column with 0
      # --> NA values occur in the sd column when length-one vector is present (genotype class count == 1)
     df.mean$sd[is.na(df.mean$sd)] <- 0 # *OBS*: df.mean$sd==NA does *NOT* work!
#     snp1  snp2	mean	sd
#     0	1	1.454427004	2.00670755
#     0	2	0.024480740	0.06535981
#     1	1	0.055989779	0.06815189
#     1	2	0.003592626	0.07120687
#     1	NA	-0.055168152	0.07441316
#     ....
    ### Remove rows with NA. That is, we do not want to plot the combination {snp1=0, snp2=NA}.
    df.mean <- df.mean[complete.cases(df.mean),]
    
    #str(df)
    #str(df.mean)
    
    ### Predicting data
    grid <- expand.grid(snp1=c(0,1,2), snp2=c(0,1,2)) # make data frame with all combinations of alleles
    grid$probe <- predict(fullmod.multiplicative, newdata=grid) # add predicted value
    ## ^ note that the "newdata" must be numeric and not factors, since the lm() was fitted to numeric values of snp1 and snp2
    grid$snp1 <- as.factor(grid$snp1); grid$snp2 <- as.factor(grid$snp2)
    #str(grid)
    
    ### Creating annotation strings
    str_pval <- sprintf("p(F-test)=%.2e", inttest_multiplicative_pval)
    str_ld <- sprintf("LD=%.2f", snp_ld)
    str_reg <- sprintf("y=%.2f + %.2f*snp1 + %.2f*snp2 + %.2f*snp1:snp2", beta_mu, beta_snp1, beta_snp2, beta_int)
    str_n_usable <- sprintf("n_usable_data_points = %d", n_usable_data_points)
    str_n_missing <- sprintf("n_missing_data_points = %d", n_missing_data_points)
    str_combined <- paste(str_pval, str_ld, str_reg, str_n_usable, str_n_missing, sep="\n")
    str_title <- sprintf("%s | %s:%s", name.probe, name.snp1, name.snp2)
    
    p <- ggplot()
    ### plot points - *removing* "NA" values: na.omit()
    p <- p + geom_point(data=na.omit(df), aes(x=snp1, y=probe, color=snp2), alpha=plot.points.alpha)
    ### mean expression
    p <- p + geom_line(data=df.mean, aes(x=snp1, y=mean, group=snp2, color=snp2))
    ### predicted
    p <- p + geom_line(data=grid, aes(x=snp1, y=probe, group=snp2, color=snp2), linetype="dashed", size=1.5)
    if (add_text_annotation) {
      ### add text annotation
      p <- p + annotate(geom="text", x=3, y=max(df$probe), label=str_combined, vjust=1.25, hjust=1)
    }
    ### Setting labs
    #p <- p + labs(title=str_title, x="SNP 1 allele", y="Probe expression (PEER residual)")
    p <- p + labs(title=str_title, x=paste0(name.snp1," allele"), y="Probe expression (PEER residual)")
    ### Setting legend
    #p <- p + guides(color=guide_legend(title="SNP 2 allele"))
    p <- p + guides(color=guide_legend(title=paste0(name.snp2,"\nallele")))

    #plots[[length(plots)+1]] <- p
    #return(plots)
    
    if (add_data_heatmap) {
      #print("YOU HAVE NOT IMPLEMENTED THIS YET!")
      
      snp1.allele <- factor(snp1,levels=c("0", "1", "2", NA), exclude=NULL)
      snp2.allele <- factor(snp2,levels=c("0", "1", "2", NA), exclude=NULL)
      #print(snp1.allele)
      #print(snp2.allele)
      tab.allele <- table(data.frame(snp1.allele, snp2.allele))
      #tab.allele
      df.twolocus <- as.data.frame(tab.allele) #   snp1.allele	snp2.allele	Freq
      ### Plot
      p.data_heatmap <- ggplot(df.twolocus, aes(x=snp1.allele, y=snp2.allele)) + geom_tile(aes(fill=Freq))
      p.data_heatmap <- p.data_heatmap + geom_text(aes(label = Freq), color="white")
      #p.data_heatmap <- p.data_heatmap + labs(x="SNP 1 allele", y="SNP 2 allele")
      p.data_heatmap <- p.data_heatmap + labs(x=paste0(name.snp1," allele"), y=paste0(name.snp2," allele"))
      p.data_heatmap <- p.data_heatmap + guides(fill=FALSE)
      
      ### List of plots
      plots <- list(p, p.data_heatmap) 
    } else {
      plots <- list(p)
    }
    
    if (save_images & (inttest_multiplicative_pval <= save_significant_treshold)) {
      print("saving plot")
      #folder_name <- "./" # OBS: remember trailing backslash "/"
      file_name <- sprintf("%d_%s|%s-%s.pdf", i, name.probe, name.snp1, name.snp2) # CONSIDER USING .SVG FORMAT
      file_name <- paste(path.out, file_name, sep="")
      
      ### Save using print()
      pdf(file_name, width=10, height=4)
      multiplot(plotlist=plots, cols=2)
      #multiplot(p.data_heatmap, cols=2)
      dev.off()
      
      ### Save using ggsave()
      #ggsave(filename=file_name, plot=p, width=10, height=6)
    }
  }
}






# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
