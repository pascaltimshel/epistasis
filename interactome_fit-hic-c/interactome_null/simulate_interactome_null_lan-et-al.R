############### SYNOPSIS ###################
## COPY OF "simulate_interactome_null.R"
## This script if meant to process the interaction table from Lan et al

### NOTES
# See Evernote for a description of the criteria to the null.

### OUTPUT
# EXPORT df.perm.concat:  a tab seperated file with indexes for the "partnerB" for the interaction.
#                         DOWNSTREAM description: interaction pairs for EXPERIMENTS are formed by joining the "partnerA" (index "a") column with the index in the null table (index "b") for the given experiment. Note that for the experiment hic_1 (positives), we always have that index a==b; whereas for null experiments we have a != b
#                             The original interactions that the NULL INDEX POINT TO (link/refer to) are found in the file, e.g.: "interation_table.fit-hi-c.nosex.interchromosomal.hESC.q_1e-13.txt"
#                         columns = experiments (hic_1, null_1, null_2,...)
#                         rows = interaction_no
#                         *REMARKS #1*: notice that the first column of df.perm.concat is the hic_1 experiment. The indexes in this column are 1, 2, 3, .., n_interactions
#                         *REMARKS #2*: remember that all experiments for all interaction pairs share the same "partnerA". Only "partnerB" is different.

###########################################

library(ggplot2)
library(permute)
library(tools) # for file_path_sans_ext()
library(dplyr)

rm(list=ls())

wd <- path.expand('~/git/epistasis/interactome_fit-hic-c/interactome_null')
setwd(wd)

###################### SOURCE ####################
source(file="function_perm_restricted.R")
##################################################


###################### PARAMS ####################
p.q.threshold <- c("OUTLIER_RM")
hic_cell_type <- "lan-et-al_K562"
n_perm <- 1000

#################################### SCRIPT #################################

######### "FOR LOOP" - START #########
time_start <- proc.time()

str.path <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/interation_table.fit-hi-c.nosex.interchromosomal.%s.q_%s.txt" # e.g. interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-07.txt
file.interaction_table <- sprintf(str.path, hic_cell_type, p.q.threshold)
file.interaction_table

df.interaction_table <- read.table(file=file.interaction_table, h=T, sep="\t", stringsAsFactors=F)

### Safety check: make sure the correct SORTED table is loaded ###
stopifnot(all(df.interaction_table$chrA <= df.interaction_table$chrB)) # --> *MUST* be TRUE!


#################### df.null ####################
#df.null <- subset(df.interaction_table, select=c("chrA", "posA", "chrB", "posB"))

#################### Run function ####################


source(file="function_perm_restricted_v2.R")
r.list <- perm_restricted(n_perm, df.interaction_table)
df.perm <- r.list[[1]]
df.perm.stats <- r.list[[2]]
######################################################


################################################################################################
###################################### Make diagnostics ########################################
################################################################################################

############ Check that Hi-C and null does not overlay INDEXES
df.tmp.check <- cbind(hic_1=seq(nrow(df.perm)), df.perm)
df.perm.hic.dup <- data.frame(n_duplicates=apply(df.tmp.check, 1, function(row) {sum(row[-1]==row[1])})) # comparing the first element (hi-c index) to all others (null index)
stopifnot(!any(df.perm.hic.dup$n_duplicates > 0))

############ Histogram sampling pool size
df.plot.pool_size <- df.perm.stats %>% select(interaction_identifer, n_pool.final) %>% arrange(n_pool.final)
df.plot.pool_size.melt <- melt(df.plot.pool_size, id.vars=c("interaction_identifer"))
n_sampling_without_replacement <- sum(!df.perm.stats$sampling_with_replacement)
text.title <- sprintf("Sampling pool size | sampling w/o replacement = %s/%s [%.2f %%]", n_sampling_without_replacement, nrow(df.interaction_table), n_sampling_without_replacement/nrow(df.interaction_table)*100 )
text.title
# ggplot
p <- ggplot(df.plot.pool_size.melt, aes(x=value)) + geom_histogram()
p <- p + geom_vline(x=n_perm, color="red", linetype="dashed")
p <- p + labs(title=text.title, x="Pool size")
p
# save
basename.plot <- sprintf("diagnostics_null_samping_pool_size_%s_nperm_%s_q_%s.pdf", hic_cell_type, n_perm, p.q.threshold)
file.plot <- file.path(dirname(file.interaction_table), basename.plot)
file.plot
ggsave( file=file.plot )


############ Histogram shrinkage of sampling pool size: facet plot
df.plot.pool_shrinkage <- df.perm.stats %>% select(interaction_identifer, s1.shrinkage, s2.shrinkage, s3.shrinkage) %>% arrange(desc(s2.shrinkage))
df.plot.pool_shrinkage.melt <- melt(df.plot.pool_shrinkage, id.vars=c("interaction_identifer"))
# ggplot
p <- ggplot(df.plot.pool_shrinkage.melt, aes(x=value)) + geom_histogram()
p <- p + labs(title="Sampling pool shrinkage", x="Pool size shrinkage")
p <- p + facet_wrap(~variable)
p
# save
basename.plot <- sprintf("diagnostics_null_samping_pool_shrinkage_%s_nperm_%s_q_%s.pdf", hic_cell_type, n_perm, p.q.threshold)
file.plot <- file.path(dirname(file.interaction_table), basename.plot)
file.plot
ggsave( file=file.plot )

############ Histogram: number of duplicated interactions across all 'experiments'
df.perm.dup <- data.frame(n_duplicates=apply(df.perm, 1, function(row) {sum(duplicated(row))})) # finding the number of indentical null for each interaction
# ggplot
ggplot(df.perm.dup, aes(x=n_duplicates)) + geom_histogram() + labs(title="Distribution of number of duplicated interactions", x="Number of duplicated interactions", y="Frequency")
# save
basename.plot <- sprintf("diagnostics_null_number_of_duplicated_interactions_%s_nperm_%s_q_%s.pdf", hic_cell_type, n_perm, p.q.threshold)
file.plot <- file.path(dirname(file.interaction_table), basename.plot)
file.plot
ggsave( file=file.plot )

############ Histogram: colSums [DOES NOT MATTER]
colsum.df.perm <- colSums(df.perm)
qplot(colsum.df.perm)


################################################################################################
############################################ EXPORT ############################################
################################################################################################
### Concatenate data frame 
df.perm.concat <- data.frame(hic_1=seq(nrow(df.perm)), df.perm)

#####################################
### Setting pathname
file.interaction_table.basename.sans_ext <- file_path_sans_ext(basename(file.interaction_table)) # removing trailing ".txt"
file.interaction_table.basename.sans_ext
file.null_table.basename.sans_ext <- sub('interation_table', 'null_table', file.interaction_table.basename.sans_ext) # regex substitution
file.null_table.sans_ext <- file.path(dirname(file.interaction_table), file.null_table.basename.sans_ext)
file.null_table <- sprintf("%s.nperm_%s.txt", file.null_table.sans_ext, n_perm)
file.null_table

### 
write.table(df.perm.concat, file=file.null_table, col.names=T, row.names=F, quote=F, sep="\t")
#?write.table

time_elapsed = proc.time() - time_start
print(sprintf("Time elapsed, q=%s: %s", p.q.threshold, time_elapsed[3]))

######### "FOR LOOP" - END #########
