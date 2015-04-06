############### SYNOPSIS ###################
# ...
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- path.expand('~/git/epistasis/interactome_fit-hic-c')
setwd(wd)

############################# READING GENE LISTs #################################
path.datafiles <- path.expand('~/Dropbox/0_Projects/p_HiC/Ferhat_Ay_2014/interaction_snpsets/maf_5_sets')

###### Read into a list of files - PATTERN VERSION - read ALL .txt files in directory:
#files <- list.files(path = path.datafiles, pattern = "*.txt", full.names = TRUE) #full path
#names(files) <- list.files(path = path.datafiles, pattern = "*.txt") # filename
#cat(names(files), sep="\n")


### From "fastEpi_compiled_broad_scp/" - UPDATED LAST TIME: 03/27/2015
# 500_snppool_hESC_q_1e-14 # hESC_width_500_maf_5_q_1e-14_epi1_1e-8
# 500_snppool_hESC_q_1e-16 # hESC_width_500_maf_5_q_1e-16_epi1_1e-8
# 1000_snppool_hESC_q_1e-12 # hESC_width_1000_maf_5_q_1e-12_epi1_1e-10
# 2500_snppool_hESC_q_1e-13 # hESC_width_2500_maf_5_q_1e-13_epi1_1e-10
# 500_snppool_hIMR90_q_1e-06 # hIMR90_width_500_maf_5_q_1e-06_epi1_1e-8
# 500_snppool_hIMR90_q_1e-08 # hIMR90_width_500_maf_5_q_1e-08_epi1_1e-8
# 2500_snppool_hIMR90_q_1e-07 # hIMR90_width_2500_maf_5_q_1e-07_epi1_1e-8
# 50000_snppool_hIMR90_q_1e-09 # hIMR90_width_50000_maf_5_q_1e-09_epi1_1e-10

filenames2read <- c("500_snppool_hESC_q_1e-14",
                    "500_snppool_hESC_q_1e-16",
                    "1000_snppool_hESC_q_1e-12",
                    "2500_snppool_hESC_q_1e-13",
                    "500_snppool_hIMR90_q_1e-06",
                    "500_snppool_hIMR90_q_1e-08",
                    "2500_snppool_hIMR90_q_1e-07",
                    "50000_snppool_hIMR90_q_1e-09")




###### Read SPECIFIC FILES:
files <- as.list(paste(path.datafiles, filenames2read, "snp_sets/set_AB.txt", sep="/"))
names(files) <- filenames2read
files
list_of_data <- llply(files, read.csv)#row.names = 1 --> NO!, stringsAsFactors = FALSE
names(list_of_data)


#################################### CREATING MATRIX #####################################
dflist <- list_of_data # copying dflist.new [TODO: make the below code a function]
mat.intersect <- matrix(data=NA, nrow=length(dflist), ncol=length(dflist), dimnames=list(names(dflist), names(dflist)))
for (i in 1:length(dflist)) {
  #gene_list = dflist[[i]][,1]
  for (j in i:length(dflist)) {
    #cat("i=",i,"j=",j,"\n")
    mat.intersect[i,j] <- length(intersect(dflist[[i]][,1], dflist[[j]][,1])) # Symmetric matrix
    mat.intersect[j,i] <- length(intersect(dflist[[i]][,1], dflist[[j]][,1])) # Symmetric matrix
    # TRY ODDs ratio
    #mat.intersect[i,j] <- length(intersect(dflist[[i]][,1], dflist[[j]][,1]))/length(dflist[[j]][,1]) # NORMALIZING by "column" length
  }
}

melt.mat.intersect <- melt(mat.intersect)

### Renaming gene lists
# rename_map <- c("gene_associated.txt"="Associated Genes", 
#                 "gene_nearest.txt"="Nearest Genes")
# levels(melt.mat.intersect$Var1)
# melt.mat.intersect$Var1 <- revalue(melt.mat.intersect$Var1, rename_map)
# melt.mat.intersect$Var2 <- revalue(melt.mat.intersect$Var2, rename_map)
levels(melt.mat.intersect$Var1)


#################################### PLOTTING #####################################

q1 <- ggplot(melt.mat.intersect, aes(Var1, Var2, fill = value)) + geom_tile()
base_size <- 12
q1 <- q1 + theme_grey(base_size = base_size) + labs(x = "", y = "", title="SNPset intersections") 
q1 <- q1 + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
#q1 <- q1 + theme(axis.text.x=element_text(size = base_size*0.8, angle=330, hjust=0, colour="grey50"))
q1 <- q1 + theme(axis.text.x=element_text(size = base_size*0.8, angle=315, hjust=0, colour="grey50"))
q1 <- q1 + geom_text(aes(label = value), size=4)
#q1 <- q1 + coord_fixed() # equal axis
q1 + scale_fill_gradient(low = "#fee8c8",  high = "#e34a33") # light_orange-->orange
#q1 + scale_fill_gradient(low = "#ece2f0",  high = "#1c9099") # baise-->green
#q1 + scale_fill_gradient(low = "#f0f0f0",  high = "#636363") # black-->white
#q1 + scale_fill_gradient(low = "blue",  high = "green")


#ggsave("SNPset_intersection-8x6.pdf", w=8, h=6)


#q1 + scale_color_brewer()
# http://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/


