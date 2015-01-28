############### SYNOPSIS ###################
# CLASS: function
# PURPOSE: write out reduced/subset of plink file
# INPUT: 
  # 1) df:        a data frame with the COLUMNS: snp1 and snp2 [more columns may exist]
  # 2) path.out:  a path to write the plink file set
# OUTPUT:
  # 1) df.snp_
  # The function prints the number of SNPs in the geno file.

# This function will identify the union of all the SNPs from the input data frame and make a new set of plink files (.bed, .bim, .fam)
# The output *directory* is path.out.
# The output *filename* is egcut_snpsubset.{bed,bim,fam,log}

### STEPS ###
# 1) Calculate SNP union
# 2) Write out "union" file. This will be the input file to PLINK
# 3) Make system call to Plink2 with the --make-bed option.
############################################


write_snp_subset_plink_file <- function(df, path.out) {
  SNP_A <- df$snp1
  SNP_B <- df$snp2
  
  
  file.out.union <- paste0(path.out,"/SNP_AB.union.txt")
  file.out.plink_prefix <-paste0(path.out,"/egcut_snpsubset")
  
  exec.plink2 <- "/Users/pascaltimshel/p_bioinformatic_tools/plink_1.9_mac_25_Nov_2014/plink"
  bed.file <- "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/all_clean/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean" # NO EXTENSION (no .bed)

  
  ### UNIQUE SNPs and probes
  SNP_A.unique <- unique(SNP_A)
  SNP_B.unique <- unique(SNP_B)  
  SNP_AB.intersect <- intersect(SNP_A.unique, SNP_B.unique)
  SNP_AB.union <- union(SNP_A.unique, SNP_B.unique)
  
  #probes.unique <- unique(df$PHENOTYPE); length(probes.unique) # CONSIDER INCLUDING THIS
  
  ### Creating data frame
  df.snpstats <- data.frame(
    n_SNP_A.unique=length(SNP_A.unique),
    n_SNP_B.unique=length(SNP_B.unique),
    n_SNP_AB.intersect=length(SNP_AB.intersect),
    n_SNP_AB.union=length(SNP_AB.union)
    )
  
  ### write UNION snps to file
  print(sprintf("Writing SNP union file: %s", file.out.union))
  cat(SNP_AB.union, sep="\n", file=file.out.union)
  
  
  ########### MAKING PLINK CALL ###########

  cmd <- sprintf("%s --bfile %s --extract %s --make-bed --out %s", exec.plink2, bed.file, file.out.union, file.out.plink_prefix)
  print(sprintf("Making plink2 cmd: %s", cmd))
  system(cmd)
  
  #system("ls -l")
  #system2("ls", "-l")
  
  return(df.snpstats)
}



