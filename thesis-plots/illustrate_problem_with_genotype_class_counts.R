############### SYNOPSIS ###################
# CLASS: script; thesis illustration
# PURPOSE: illustrate the problem with low genotype class counts

# *IMPORTANT*: this script is a PARTIAL COPY of another script "explore_epistasis_table.R"

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


################################################################################################
###################################### Settings/params ########################################
################################################################################################

### Job
#param.job_identifier <- "hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8" # n=38,987
#param.job_identifier <- "hIMR90_width_1000_maf_5_q_1e-06_epi1_1e-10" # n=96,083 | EIID_probe
param.job_identifier <- "hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8" # n=102,339 | EIID_probe

#param.job_identifier <- "hESC_width_2500_maf_5_q_1e-13_epi1_1e-10" # n=248,510

### Epistasis table
# this is needed for loading the correct epistasis table
#param.epistasis_table.type <- "epistasis_table.txt"
#param.epistasis_table.type <- "epistasis_table_pruned_hemani.txt"
param.epistasis_table.type <- "epistasis_table_pruned_EIID_probe.txt"

### Paths - out
path.out.subdir <- "epistasis_table_processing"

path.out <- path.expand("~/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2")
path.out.base <- file.path(path.out, param.job_identifier, path.out.subdir)
path.out.base.analysis_epi_table <- file.path(path.out.base, "analysis_epi_table_GENOTYPE_CLASS_COUNT") #***OBS*** NEW
path.out.base.analysis_epi_table

path.out.base.plots_model <- file.path(path.out.base.analysis_epi_table, "plot_epistasis_model")

path.out.base.PLINK_files <- path.out.base.analysis_epi_table

### Paths - libary (scripts)
path.library <- path.expand("~/git/epistasis/interactome_fit-hic-c/explore_epi_model")

### Files
file.epi_table <- file.path(path.out.base, param.epistasis_table.type)
file.epi_table # e.g. "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/fastEpi_compiled_broad_scp_null_v2/hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8/epistasis_table_pruned_hemani.txt"



################################################################################################
######################################## Initialize ##########################################
################################################################################################

### Source ###
source(file.path(path.library,"function_write_snp_subset_plink_file.R"))
# write_snp_subset_plink_file()
source(file.path(path.library,"function_min_genotype_class_count.R"))
# subset_genotype_class_count()
source(file.path(path.library,"function_Replication_fastEpistasis_comparison.R"))
# RunEpistasisTests()
source(file.path(path.library,"function_plot_epistasis_model.R"))
# plot_EpiModel()


### Create directories if needed
if (file.exists(path.out.base.analysis_epi_table)) {
  print("OBS: path.out.base.analysis_epi_table folder exists")
}
dir.create(path.out.base.analysis_epi_table, showWarnings=FALSE) # dir.create() does NOT crash if the directory already exists, it just prints out a warning.
print(path.out.base.analysis_epi_table) # no trailing slash


################################# LOAD epistasis_table ############################
### EXAMPLE - April 2015 ###
### REMARKS
# - no header
# 5  rs7722032  21  rs7281390	-0.23753	57.45355	3.46056E-14	ILMN_1699160	hic_1_5965	True
# 1	rs11122320	2	rs4954710	+0.12256	53.48095	2.61113E-13	ILMN_1688234	hic_1_2276	True
# ...
# 4  rs7676504	18	rs1943313	+0.33703	58.93330	1.63108E-14	ILMN_1791759	null_1_7165	True
# 12	rs7966152	11	rs968064	+0.32100	53.49671	2.59026E-13	ILMN_2065022	null_1_2091	True

### READ table
#epi_table.header <- c("CHR","SNP_A","CHR","SNP_B","BETA","CHISQ","PVALUE","PHENOTYPE","EIID","SIGNIFICANCE") # KEEP THIS HEADER UPDATED!
#epi_table.header <- c("chr1", "snp1", "chr2", "snp2", "BETA", "CHISQ", "PVALUE", "probename","EIID","EID",) # KEEP THIS HEADER UPDATED!
epi_table.header <- c("chr1", "snp1", "chr2", "snp2", "BETA", "CHISQ", "PVALUE", "probename","EIID","EID",
                      "RefSeq_ID","Symbol","Chromosome","Probe_Chr_Orientation","Probe_Coordinates_start","Probe_Coordinates") # KEEP THIS HEADER UPDATED!
df.epi_table <- read.delim(file.epi_table, stringsAsFactors=F, header=F, col.names=epi_table.header)#
str(df.epi_table)
nrow(df.epi_table)

### ***Creating EID column*** 
#[*OUTCOMMENT ME LATER*]
# df.epi_table$EID <- sub("_\\d*?$", "", df.epi_table$EIID, perl=TRUE) # OBS: meta characters need to be escaped (e.g. \d --> \\d)
## REF for R regex: http://www.johndcook.com/blog/r_language_regex/


### UNIQUE SNPs
SNP_AB.union <- union(df.epi_table$snp1, df.epi_table$snp2); length(SNP_AB.union)


######################################## *CREATE df.enrichment* ###########################################
EIDs <- c("hic_1", paste0("null_",seq(1,1000)))
df.enrichment <- data.frame(EID=EIDs)
#df.enrichment
### Calculate step enrichment
df.enrichment.1.no_min_geno <- df.epi_table %>% group_by(EID) %>% dplyr::summarise(count=n()) %>% arrange(desc(count))
df.enrichment.1.no_min_geno <- dplyr::rename(df.enrichment.1.no_min_geno, count1_no_min_geno = count) # rename variable
df.enrichment.1.no_min_geno
### UPDATE df.enrichment
df.enrichment <- full_join(df.enrichment, df.enrichment.1.no_min_geno, by="EID")



######################################## Read probes ###########################################
### Read file with 9269 probes and 832 individuals
#file.probes <- "/Users/pascaltimshel/Dropbox/Ubuntu/peer_analysis/export_residuals/egcut.peer_residuals_log2_k50.top50_mean_top50_var_refseq.txt"
# df.probes <- read.table(file.probes, h=T)
# rownames(df.probes) <- df.probes$IID
# df.probes <- subset(df.probes, select=-c(IID, FID))
# save(df.probes, file="RData_tmp/tmp_EGCUT_9269_probes.RData")
load(file="/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/R_epistasis_tests/RData_tmp/tmp_EGCUT_9269_probes.RData") # df.probes

###############################################################################################################
######################################## Read FastEpistasis results ###########################################
###############################################################################################################

### *BACKWARDS COMPATIBILITY* ###
# We are using the data frame with variable name "df.fastEpistasis.results" because of legacy of previously written code
df.fastEpistasis.results <- df.epi_table

### Create ID column: probe.snp1.snp2
df.fastEpistasis.results$ID <- apply(df.fastEpistasis.results[,c("probename", "snp1", "snp2")], 1, function(df.row) {paste(df.row[1], paste(sort(c(as.character(df.row[2]), as.character(df.row[3]))),collapse="."),sep=".")})
str(df.fastEpistasis.results)

###############################################################################################################
######################################## Write PLINK files ###########################################
###############################################################################################################

### WRITING NEW SET OF PLINK FILES (from function_write_snp_subset_plink_file.R)
# --> output set of plink files will be written to PATH.OUT/egcut_snpsubset.{bed,bim,fam,log}
df.snpstats <- write_snp_subset_plink_file(df=df.fastEpistasis.results, path.out=path.out.base.PLINK_files)
df.snpstats

######################################## Read Genotypes ###########################################
#### Small test set
#bed.file <- "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/SNPs12959/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.fit-hi-c.SNPs12959"
#bed.file <- "/Users/pascaltimshel/Dropbox/5_Data/EGCUT_DATA/geno/all_clean/Prote_370k_251011.no_mixup.with_ETypes.chr_infered.clean"
bed.file <- paste0(path.out.base.PLINK_files, "/egcut_snpsubset")
#### Read data:
sample <- read.plink(bed.file)
# Converting SnpMatrix to Numeric Matrix
geno <- as(sample$genotypes, "numeric") # 0,1,2. Missing data will be NA.
dim(geno)

######################################## File checks ###########################################
df.check.sync.probe_geno <- data.frame(rownames(df.probes), rownames(geno))
all(df.check.sync.probe_geno[,1]==df.check.sync.probe_geno[,2]) # --> TRUE: in sync



######################################## Calculate MAF ###########################################
### snpStats - MAF for all SNPs in the loaded genotype file
snpsum <- col.summary(sample$genotypes)
df.snpsum.maf <- data.frame(snp=rownames(snpsum), MAF=snpsum$MAF)
head(df.snpsum.maf)

### Adding column to data frame: MAF
df.fastEpistasis.results$snp1.maf <- df.snpsum.maf$MAF[match(df.fastEpistasis.results$snp1, df.snpsum.maf$snp)]
df.fastEpistasis.results$snp2.maf <- df.snpsum.maf$MAF[match(df.fastEpistasis.results$snp2, df.snpsum.maf$snp)]
### Alternative approach: LOOP OVER VECTOR (sapply()) - MUCH SLOWER!
#snp1.maf <- sapply(df.fastEpistasis.results$snp1, function(rsID) {df.snpsum.maf$MAF[match(rsID, df.snpsum.maf$snp)]})
str(df.fastEpistasis.results)

### Adding column to data frame: pairwise MAF
df.fastEpistasis.results$snp.maf.min <- with(df.fastEpistasis.results, pmin(snp1.maf, snp2.maf))
#sum(df.fastEpistasis.results$snp.maf.min > 0.2)

df.fastEpistasis.results.maf.gt <- subset(df.fastEpistasis.results, snp.maf.min>=0.04) # *OBS: only TEMPORARY*
#df.fastEpistasis.results.maf.gt <- subset(df.fastEpistasis.results, snp.maf.min>=0.1)
#df.fastEpistasis.results.maf.gt <- subset(df.fastEpistasis.results, snp.maf.min>=0.25)
nrow(df.fastEpistasis.results.maf.gt)

#################################### Calculate min class count ######################################

source(file.path(path.library,"function_min_genotype_class_count.R"))
# From script: function_subset_min_genotype_class_count.R
df.input <- df.fastEpistasis.results.maf.gt
#df.input <- df.fastEpistasis.results
#idx <- 1:1000
idx <- 1:nrow(df.input)
df.min_genotype_class_count.full <- min_genotype_class_count(idx=idx, df.input=df.input)
nrow(df.min_genotype_class_count.full)
# hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8 | min 3 --> n=703

### *THRESHOLD* ###
min_genotype_class_count_threshold <- 0 # 3 | 5 | 10

### SUBSET DATA FRAME
df.min_genotype_class_count.sub <- subset(df.min_genotype_class_count.full, min_genotype_class_count >= min_genotype_class_count_threshold)
print(sprintf("number of SNP-pairs satisfying criteria: %s (%.2f %%)", nrow(df.min_genotype_class_count.sub), nrow(df.min_genotype_class_count.sub)/nrow(df.min_genotype_class_count.full)*100))

########## TEMPORARY - SAVE ########
#save(df.min_genotype_class_count.full, file="df.min_genotype_class_count.full.hIMR90_l-2500.RData")
#save(df.min_genotype_class_count.full, file="df.min_genotype_class_count.full.hIMR90_l-1000.RData")

################################ *** EPISTASIS ENRICHMENT *** #################################

#############################
### Calculate step enrichment
df.enrichment.2.min_geno <- df.min_genotype_class_count.sub %>% group_by(EID) %>% dplyr::summarise(count=n()) %>% arrange(desc(count))
df.enrichment.2.min_geno <- dplyr::rename(df.enrichment.2.min_geno, count2_min_geno = count) # rename variable
df.enrichment.2.min_geno
### UPDATE df.enrichment
df.enrichment <- full_join(df.enrichment, df.enrichment.2.min_geno, by="EID")

### SUBSTITUTE ALL "NA" with 0
df.enrichment[is.na(df.enrichment)] <- 0
#############################

########### ***CALCULATE P-value*** ###############
df.tmp.enrichment.hic <- df.enrichment %>% filter(EID=="hic_1")
df.tmp.enrichment.hic
sum(df.enrichment$count1_no_min_geno >= df.tmp.enrichment.hic$count1_no_min_geno)/1000
sum(df.enrichment$count2_min_geno >= df.tmp.enrichment.hic$count2_min_geno)/1000


##################################### *** ARRANGE/SORT the EPISTASIS table ***  ###########################################
str(df.min_genotype_class_count.sub)
df.x <- df.min_genotype_class_count.sub %>% arrange(desc(min_genotype_class_count),desc(BETA))
#df.x <- df.min_genotype_class_count.sub %>% arrange(desc(BETA), desc(min_genotype_class_count))
nrow(df.x)

#####################################################################################################
##################################### Plot epistasis model  #########################################


#ILMN_2065022 |rs10041179, rs9294148 #--> USE THIS | min_geno = 0 | PVALUE 2e-110
#plot_epistasis_model3_ILMN_1680436|rs1478590-rs10132344 #--> USE THIS | min_geno = 1
#plot_epistasis_model8_ILMN_1781285|rs322673-rs2329119 #--> USE THIS | min_geno = 9

#plot_epistasis_model46_ILMN_2385239|rs4519460-rs1224449 # min_geno = 6
df.selected <- df.min_genotype_class_count.sub %>% filter(
                                          (snp1=="rs10041179" & snp2=="rs9294148" & probename=="ILMN_2065022")
                                          |(snp1=="rs1478590" & snp2=="rs10132344" & probename=="ILMN_1680436") 
                                          |(snp1=="rs322673" & snp2=="rs2329119" & probename=="ILMN_1781285")
                                          )
df.selected


idx <- 1:10
#idx <- 1:nrow(df.probe_snp_pair)
### Run function from script "function_plot_epistasis_model.R"
source(file.path(path.library,"function_plot_epistasis_model.R"))
#plots <- plot_EpiModel(idx, df.input=df.x, path.out=path.out.base.plots_model, add_data_heatmap=TRUE, save_images=TRUE, save_significant_treshold=1e-5, plot.points.alpha=0.7)
plots <- plot_EpiModel(1:nrow(df.selected), df.input=df.selected, path.out=path.out.base.plots_model, add_text_annotation=FALSE, add_data_heatmap=TRUE, save_images=TRUE, save_significant_treshold=1e-5, plot.points.alpha=0.7)

###########################################################################################
############################# PLOT P-value vs MIN-MAF/GENO ################################
###########################################################################################

df.tmp.y <- df.min_genotype_class_count.full %>% arrange(PVALUE) %>% slice(1:5)
df.tmp.y

### Draw random SNP-probe pairs
set.seed(1) # important for reproduction
df.rand <- df.min_genotype_class_count.full %>% sample_n(100000)
#df.rand <- df.min_genotype_class_count.full
str(df.rand)

with(df.rand, cor(min_genotype_class_count, -log10(PVALUE)))
with(df.rand, cor(min_genotype_class_count, PVALUE))
with(df.rand, cor(snp.maf.min, PVALUE))

### correlation
ggplot(df.rand) + geom_point(aes(x=snp.maf.min, y=min_genotype_class_count))
with(df.rand, cor(min_genotype_class_count, snp.maf.min))
with(df.rand, cor.test(min_genotype_class_count, snp.maf.min, method="pearson"))

### MIN maf
ggplot(df.rand) + geom_point(aes(x=snp.maf.min, y=PVALUE))
ggplot(df.rand) + geom_point(aes(x=snp.maf.min, y=-log10(PVALUE)))

### MIN genotype class count
ggplot(df.rand) + geom_point(aes(x=min_genotype_class_count, y=-log10(PVALUE)))


#### USE this one for THESIS ####
p <- ggplot(df.rand)
p <- p + geom_boxplot(aes(x=factor(min_genotype_class_count), y=-log10(PVALUE)))
p <- p + labs(x="Minimum genotype class count", y=expression(paste(-log[10], ("p-value") ) ))
p
ggsave(file="minimum_genotype_class_count-vs-pvalue_boxplot-8x4.pdf", w=8, h=4)




###########################################################################################
############################# Barplot | distribution of minimum genotype class count ################################
###########################################################################################

### ** FINISH MEEE **** Stopped here July 10th ###

#df.table <- data.frame(table(df.rand$min_genotype_class_count))
#str(df.table)
#df.table %>% rename(Var1)

df.table <- df.rand %>%
  group_by(min_genotype_class_count) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n), percentage = freq*100)
df.table

p <- ggplot(df.table)
p <- p + geom_bar(aes(x=min_genotype_class_count, y=freq), stat="identity")
p <- p + labs(x="Minimum genotype class count", y="Frequency")
p + scale_x_continuous(breaks=0:6, limits=c(-0.5,6))
p

ggsave(file="minimum_genotype_class_count_counts_barplot_xlim-8x4.pdf", w=8, h=4)

###########################################################################################




