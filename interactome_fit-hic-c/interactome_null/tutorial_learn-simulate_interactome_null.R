#install.packages("gtools")
#install.packages("permute")
#install.packages("dplyr")

library(ggplot2)
library(permute)
library(gtools)
library(plyr)
#library(dplyr)
#?dplyr


rm(list=ls())

set.seed(1)

x <- 1:100
n_perm <- 1000 # 100000
mat <- t(replicate(n_perm, sample(x)))

### Calculate "error rate" for each position
x.err <- sapply(seq_along(x), function(i) {sum(mat[,i]==x[i])/dim(mat)[1]})
#x.err1 <- numeric(); for (i in seq_along(x)) {x.err1[i] <- sum(mat[,i]==x[i])} # ALTERNATIVE SOLUTION
x.err # --> the error rate is 1/x
sprintf("Mean Error Rate: %s", mean(x.err))
sprintf("1/x: %s", 1/length(x))

df <- rdply(n_perm, sample(x)) # first column in this matrix (df$.n) is the original index



###################### USING "permute" library #################
#shuffle(20) # --> generates vector
?shuffleSet
df.shuffle <- shuffleSet(100, 5) # shuffleSet(n, nset)
#class(df.shuffle) # --> "permutationMatrix" "matrix"
df.shuffle <- as.data.frame(t(df.shuffle))

?how
CTRL <- how(within = Within(type = "series", mirror = TRUE))
#CTRL <- how(within = Within(type = "grid", mirror = TRUE)) # --> does not work: Set of permutations < 'minperm'. Generating entire set.
#CTRL <- how(within = Within(type = "free", mirror = TRUE)) # -->
CTRL
shuffleSet(5, 10, CTRL)

#Permute package
#http://www.r-bloggers.com/generating-sets-of-permutations/
#shuffle(20)
#shuffleSet(3, 10)


#How to resample in R without repeating permutations?
#http://stats.stackexchange.com/questions/24300/how-to-resample-in-r-without-repeating-permutations

#Generating all distinct permutations of a list in R
#http://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r

is.language(call("sample", 1))
call("sample", 1)
is.call()
#http://www.pmc.ucsc.edu/~mclapham/Rtips/resampling.htm

#If expr is a function call, be aware of assumptions about where it is evaluated, and in particular what ... might refer to. You can pass additional named arguments to a function call as additional named arguments to replicate: see "Examples".