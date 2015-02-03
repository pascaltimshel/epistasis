############### SYNOPSIS ###################
# CLASS: script; data/table generative
# PURPOSE: Generate NULL interactions with as few "true positives" as possible

# STEPS
# 1) define function with restrictions on permutations
# 2) run function
# 3) (make a few plots)

###########################################

library(ggplot2)
library(permute)
library(tools) # for file_path_sans_ext()

rm(list=ls())

wd <- '/cvar/jhlab/timshel/git/epistasis/interactome_fit-hic-c/interactome_null'
set(wd)

###################### TODO ####################
#1) *SEMI IMPORTANT*: fix the script to take into account that some positions on the genome participate in several interactions.
  # --> this means that the current implementation may contain a few "true positives".
  # --> to fix this, load in "chr:pos" in the table and permute these instead.
  # COMPLETED this task 01/28/2015
#2) find out how many theoretical "restricted permutations", n_perm_max, that exists as a function of length(x). 
  # --> Print warning if n_perm > n_perm_max

###################### Read interaction table ####################
str.path <- "/Users/pascaltimshel/p_HiC/Ferhat_Ay_2014/interaction_tables/interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_%s.txt" # e.g. interation_table.fit-hi-c.nosex.interchromosomal.hIMR90.q_1e-07.txt
p.q.threshold <- 1e-06
#p.q.threshold <- 1e-07
#p.q.threshold <- 1e-08
#p.q.threshold <- 1e-09
#p.q.threshold <- 1e-10
file.interaction_table <- sprintf(str.path, p.q.threshold)
file.interaction_table

df.interaction_table <- read.table(file=file.interaction_table, h=T, sep="\t", stringsAsFactors=F)

### Safety check: make sure the correct SORTED table is loaded ###
all(df.interaction_table$chrA <= df.interaction_table$chrB) # --> *MUST* be TRUE!



###################### USING own function - WITH restrictions #################

perm_restricted <- function(x, n_perm) {
  ### DESCRIPTION
  ### Input
  #  x: vector to sample from. The values in x will be used to avoid any "positives".
  
  ### For debugging - running for-loop manually
  #n_perm <- 10
  
  ### Initialyzing data frames
  df.perm <- data.frame(matrix(NA, ncol=n_perm, nrow=length(x)))
  colnames(df.perm) <- paste0("null_", seq(n_perm))
  ## COLUMN=permutation_no (permuted sample)
  ## ROW="observation"
  df.perm.stats <- data.frame(attempts=rep(NA, n_perm))
  
  for (i in 1:n_perm) {
    n_perm_attempts <- 1
    flag_make_new_permutation <- TRUE
    while (flag_make_new_permutation) {
      perm <- sample.int(length(x)) # this will sample integers in the range 1...n_perm.
      #perm <- sample(x) # this will sample the values in x
      perm_val <- x[perm] # now reorder x according to the permuted index in perm.
      
      bool_identical_positions <- x == perm_val
      print(sprintf("i=%s | Attempt=%s | %s identical positions in x and perm", i, n_perm_attempts, sum(bool_identical_positions)))
      flag_make_new_permutation <- any(bool_identical_positions) # value will be TRUE if any positions are identical    
      n_perm_attempts <- n_perm_attempts + 1
    }
    print(sprintf("i=%s | Ok, found perm", i))
    df.perm[,i] <- perm
    df.perm.stats[i,] <- n_perm_attempts
  }
  df.perm.stats$error <- sapply(seq_along(df.perm), function(i) {sum(df.perm[,i]==x)/dim(df.perm)[1]})
  
  r.list <- list(df.perm, df.perm.stats)
  return(r.list)
}

set.seed(1) # Important to set seed for REPRODUCIBLE results
#x <- 1:1000
#n_perm <- 2000 # 100000
x <- paste0(df.interaction_table$chrB,":", df.interaction_table$posB)
#n_perm <- 10
n_perm <- 1000


### Run function
r.list <- perm_restricted(x, n_perm)
df.perm <- r.list[[1]]
df.perm.stats <- r.list[[2]]

###### Make diagnostics ######
### Attempts histogram
qplot(df.perm.stats$attempts, geom='histogram') + labs(title="Distribution of number of attempts")

### Number of duplicated interactions across all 'experiments'
df.perm.dup <- data.frame(n_duplicates=apply(df.perm, 1, function(row) {sum(duplicated(row))})) # finding the number of indentical null for each interaction
#x <- df.perm.dup[df.perm.dup>0]
ggplot(df.perm.dup, aes(x=n_duplicates)) + stat_bin(aes(y=..count../sum(..count..))) + labs(title="Distribution of number of duplicated interactions", x="Number of duplicated interactions", y="Frequency")


###################### EXPORT ######################
### Concatenate data frame 
df.perm.concat <- data.frame(hic_1=seq(x), df.perm)

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

