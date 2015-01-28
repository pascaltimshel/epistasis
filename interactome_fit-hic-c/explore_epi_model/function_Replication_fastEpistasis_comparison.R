############### SYNOPSIS ###################
# *IMPORTANT*: this script is a MODIFIED DUPLICATE of the "..XXX/EGCUT_DATA/R_epistasis_tests" directory


# CLASS: function
# PURPOSE: 
# INPUT: data frame with the COLUMNS: snp1, snp2, probename [more columns may exist]
# OUTPUT: a NEW data frame with columns like F_stat, P-val and more.
  # This data frame will also contain the allele count data for each SNP.
# USAGE: the main function "RunEpistasisTests" needs the expression data, genotype data and a data frame of SNPs to loop over
############################################


RunEpistasisTests <- function(idx, df.probe_snp_pair, geno, expression) { 
  # idx: e.g. 1:10
  # df.probe_snp_pair is a data frame with the columns: snp1, snp2, probename 
  
  df.res <- data.frame()
  for(i in idx) {
    cat(i, "of", length(idx), "\n")
    name.snp1 <- df.probe_snp_pair$snp1[i]
    name.snp2 <- df.probe_snp_pair$snp2[i]
    name.probe <- df.probe_snp_pair$probename[i]
    
    df.res.tmp <- EpistasisTests(name.snp1, name.snp2, name.probe, geno, expression)
    df.res <- rbind(df.res, df.res.tmp) # consider using plyr (rbind.fill) to allow for unequal number of columns
  }
  return(df.res)
}

EpistasisTests <- function(name.snp1, name.snp2, name.probe, geno, expression) {

  # Extract data
  snp1 <- geno[, colnames(geno) == name.snp1]
  snp2 <- geno[, colnames(geno) == name.snp2]
  probe <- expression[, colnames(expression) == name.probe]
  
  if ( (length(snp1)<1) | (length(snp2)<1) | (length(probe)<1) ) {
    #print(paste("length(snp1):", length(snp1)))
    #print(paste("length(snp2):", length(snp2)))
    #print(paste("length(probe):", length(probe)))
    #stop("snp1, snp2 or probe not found")
    print("snp1, snp2 or probe not found")
    return(data.frame()) # returning empty data frame
  } 

  n_usable_data_points <- sum(!is.na(snp1) & !is.na(snp2))
  n_missing_data_points <- sum(is.na(snp1)) + sum(is.na(snp2))
  
  
  ### MAF calculation
  # REMEMBER to remove NA when using sum()
  # number_of_individuals_with_snp/number_of_non_missing_datapoints
  # number_of_SNP[x]_alleles/total_number_of_alleles
  tmp.snp1.maf <- (sum(snp1==1,na.rm=T)*1 + sum(snp1==2, na.rm=T)*2)/(2*sum(!is.na(snp1)))
  tmp.snp2.maf <- (sum(snp2==1,na.rm=T)*1 + sum(snp2==2, na.rm=T)*2)/(2*sum(!is.na(snp2)))
  snp1.maf <- ifelse(tmp.snp1.maf <= 0.5, tmp.snp1.maf, 1-tmp.snp1.maf)
  snp2.maf <- ifelse(tmp.snp2.maf <= 0.5, tmp.snp2.maf, 1-tmp.snp2.maf)
  
  snp1.allele <- factor(snp1,levels=c("0", "1", "2", NA), exclude=NULL)
  snp2.allele <- factor(snp2,levels=c("0", "1", "2", NA), exclude=NULL)
  #print(snp1.allele)
  #print(snp2.allele)
  tab.allele <- table(data.frame(snp1.allele, snp2.allele))
  tab.allele
  df.twolocus <- as.data.frame(tab.allele)
  rownames(df.twolocus) <- paste("count_SNP1/SNP2=", df.twolocus$snp1.allele, "/", df.twolocus$snp2.allele, sep="")
  df.twolocus <- subset(df.twolocus, select=c(-snp1.allele, -snp2.allele))
  df.twolocus.t <- t(df.twolocus)

  
  ################## Statistical tests #####################
  ##### Multiplicative model #####
  ## Full
  fullmod.multiplicative <- lm(probe ~ snp1 + snp2 + snp1:snp2)
  fullmod.multiplicative
  summary(fullmod.multiplicative)
  beta_int_multiplicative <- fullmod.multiplicative$coefficients[4] # coefficient for interaction term (snp1:snp2)
  ## Marginal
  margmod.multiplicative <- lm(probe ~ snp1 + snp2)
  ## Test for model reduction
  inttest_multiplicative <- anova(margmod.multiplicative, fullmod.multiplicative)
  
  ### Extractions
  F_statistic <- inttest_multiplicative$F[2]
  F_df_denominator = margmod.multiplicative$df.residual - fullmod.multiplicative$df.residual # m_full - m_red = N_param_full - N_param_red = (3+1)-(2+1)=1 <-- This number will always be the same
  F_df_numerator = fullmod.multiplicative$df.residual # n - m_full = N_observations - N_param_full = 6-(3+1)=2 <--- N_observations will be *VARIABLE*
  # F(F_statistic, df_denominator, df_numerator, lower.tail=F) --> P-values
  
  #F_pval <- -log10(inttest_multiplicative$P[2]) 
  F_pval <- inttest_multiplicative$P[2] 
  
  ################## Make data frame #####################
  df.res.tmp <- data.frame(
    snp1=name.snp1, snp2=name.snp2, probename=name.probe,
    beta=beta_int_multiplicative,
    F_pval,
    F_statistic,
    F_df_denominator,
    F_df_numerator,
    snp1_maf=snp1.maf,
    snp2_maf=snp2.maf,
    n_usable_data_points,
    n_missing_data_points
    )
  
  df.res.tmp <- cbind(df.res.tmp, df.twolocus.t)
  
  return(df.res.tmp)
}


#http://stackoverflow.com/questions/18790458/how-to-attach-variables-of-one-observation-to-another-within-the-same-dataset-in
#data <- merge(data, spouses.yearly, by=c("pID", "year"), all.x=TRUE)

# sig[["count_SNP1/SNP2=0/0"]] <- NA
# sig[["count_SNP1/SNP2=1/0"]] <- NA
# sig[["count_SNP1/SNP2=2/0"]] <- NA
# sig[["count_SNP1/SNP2=NA/0"]] <- NA
# sig[["count_SNP1/SNP2=0/1"]] <- NA
# sig[["count_SNP1/SNP2=1/1"]] <- NA
# sig[["count_SNP1/SNP2=2/1"]] <- NA
# sig[["count_SNP1/SNP2=NA/1"]] <- NA
# sig[["count_SNP1/SNP2=0/2"]] <- NA
# sig[["count_SNP1/SNP2=1/2"]] <- NA
# sig[["count_SNP1/SNP2=2/2"]] <- NA
# sig[["count_SNP1/SNP2=NA/2"]] <- NA
# sig[["count_SNP1/SNP2=0/NA"]] <- NA
# sig[["count_SNP1/SNP2=1/NA"]] <- NA
# sig[["count_SNP1/SNP2=2/NA"]] <- NA
# sig[["count_SNP1/SNP2=NA/NA"]] <- NA