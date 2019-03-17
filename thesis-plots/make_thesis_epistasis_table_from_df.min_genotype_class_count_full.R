############### SYNOPSIS ###################
# CLASS: script; thesis illustration
# PURPOSE: illustrate the problem with low genotype class counts

# *IMPORTANT*: this script is a PARTIAL COPY of another script
  # ---> "illustrate_problem_with_genotype_class_counts.R"



# 2015-11-19 REMARK: THIS IS A VERY MESSY SCRIPT. YOU MAY DELETE IT OR CLEAN IT LATER....


############################################

library(snpStats)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr) # REMEMBER 
library(tools)

rm(list=ls())

wd <- path.expand("~/git/epistasis/thesis_plots")
setwd(wd)
######################

load("df.min_genotype_class_count.full.hIMR90_l-2500.RData") # --> df.min_genotype_class_count.full

################################################################################################
######################################## Read Genotypes ##########################################
################################################################################################
#### Small test set
#bed.file <- "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/SNPs12959/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.fit-hi-c.SNPs12959"
#bed.file <- "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/all_clean/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean"
file.base <- "/Users/pascaltimshel/Dropbox/0_Projects/git/epistasis/thesis_plots/PLINK_GENOTYPES-genotype_class_count_full_hIMR90_l-2500/egcut_snpsubset"

### REFERENCE: http://www.stat-gen.org/tut/tut_intro.html
gwas.fn <- lapply(c(bed='bed',bim='bim',fam='fam'), function(n) sprintf("%s.%s", file.base, n))
gwas.fn
#### Read data:
geno <- read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam)
geno

#Obtain the SNP information from geno list
genoBim <- geno$map
colnames(genoBim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(genoBim))

############# JOIN ###############
?semi_join

str(df.min_genotype_class_count.full)
## ROUND 1: SNP1
df.thesis.table <- df.min_genotype_class_count.full %>% left_join(genoBim, by=c("snp1"="SNP"))
df.thesis.table$snp1.pos <- df.thesis.table$position # Copy column, because "position" will now be overwritten
str(df.thesis.table)
df.thesis.table <- df.thesis.table %>% left_join(genoBim, by=c("snp2"="SNP"))
df.thesis.table$snp2.pos <- df.thesis.table$position.y # OBS: beware of the name! Copy column

### Probe co-localization
df.thesis.table$probe_colocalization <- with(df.thesis.table, (Chromosome == chr1) | (Chromosome == chr2))

################################################################################################
###################################### TABLE 1: "REAL EPISTASIS" with GCC>5 #####################
################################################################################################

############## SUBSET min class count ##########
df.thesis.table.true_epistasis <- df.thesis.table %>% filter(min_genotype_class_count>=10)
sum(df.thesis.table.true_epistasis$probe_colocalization)/length(df.thesis.table.true_epistasis) # --> 0.05882353
write.csv(df.thesis.table.true_epistasis, file="df.thesis.table.true_epistasis_minGCC-10.csv")


################################################################################################
###################################### TABLE 2: SPATIAL EPISTASIS HITS #####################
################################################################################################

df.thesis.table.hic_1.minGCC3 <- df.thesis.table %>% filter(EID=="hic_1", min_genotype_class_count>=3)

### Rearrangeing and formatting
df.thesis.table.hic_1.minGCC3.format <- with(df.thesis.table.hic_1.minGCC3, data.frame(
  SNP1=paste0(chr1,":",snp1.pos),
  rsID.1=snp1,
  MAF1=snp1.maf,
  SNP2=paste0(chr2,":",snp2.pos),
  rsID.2=snp2,
  MAF2=snp2.maf,
  "P-value"=PVALUE,
  Beta=BETA,
  minGCC=min_genotype_class_count,
  Gene=Symbol,
  Probe=paste0(Chromosome, ":", sapply(strsplit(df.thesis.table.hic_1.minGCC3$Probe_Coordinates, "-"),'[',1)), # extracting first element from coordinate
  EIID=EIID,
  EID=EID
))

#write.csv(df.thesis.table.hic_1.minGCC3.format, file="df.thesis.table.format.hic_1.minGCC3.csv")

################################################################################################
################# TABLE 3: EXAMPLES OF "DIRTY" EPISTASIS (MATCHING FIGURE) #####################
################################################################################################

#ILMN_2065022 |rs10041179, rs9294148 #--> USE THIS | min_geno = 0 | PVALUE 2e-110
#plot_epistasis_model3_ILMN_1680436|rs1478590-rs10132344 #--> USE THIS | min_geno = 1
#plot_epistasis_model8_ILMN_1781285|rs322673-rs2329119 #--> USE THIS | min_geno = 9

#plot_epistasis_model46_ILMN_2385239|rs4519460-rs1224449 # min_geno = 6
df.thesis.table.selected <- df.thesis.table %>% filter(
  (snp1=="rs10041179" & snp2=="rs9294148" & probename=="ILMN_2065022")
  |(snp1=="rs1478590" & snp2=="rs10132344" & probename=="ILMN_1680436") 
  |(snp1=="rs322673" & snp2=="rs2329119" & probename=="ILMN_1781285")
)

df.thesis.table.selected


# df.thesis.table.selected.format_3SNP-probe-pairs-matching-example-figure
df.thesis.table.selected <- df.thesis.table %>% filter(
  (snp1=="rs10041179" & snp2=="rs9294148")
  |(snp1=="rs1478590" & snp2=="rs10132344") 
  |(snp1=="rs322673" & snp2=="rs2329119")
)
df.thesis.table.selected


### Rearrangeing and formatting
df.thesis.table.selected.format <- with(df.thesis.table.selected, data.frame(
  SNP1=paste0(chr1,":",snp1.pos),
  rsID.1=snp1,
  MAF1=snp1.maf,
  SNP2=paste0(chr2,":",snp2.pos),
  rsID.2=snp2,
  MAF2=snp2.maf,
  "P-value"=PVALUE,
  Beta=BETA,
  minGCC=min_genotype_class_count,
  Gene=Symbol,
  Probe=paste0(Chromosome, ":", sapply(strsplit(df.thesis.table.selected$Probe_Coordinates, "-"),'[',1)), # extracting first element from coordinate
  EIID=EIID,
  EID=EID
))


### WRITE
# write.csv(df.thesis.table.selected.format, "df.thesis.table.selected.format_3SNP-probe-pairs-matching-example-figure.csv")

















########################################### THE BELOW CODE MAY BE DELETED #############################
