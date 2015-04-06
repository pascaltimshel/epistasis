############### SYNOPSIS ###################
# CLASS: function
# PURPOSE: this function sorts interactions so chr(interactionA) <= chr(interactionB)
# VALIDITY: fit-hi-c interaction tables
# INPUT: the function should be parsed a data frame with the following columns in THIS ORDER: c("chrA", "posA", "chrB", "posB")
  # Chr columns MUST be numeric!
# OUTPUT: a sorted data frame.
############################################
library(plyr)

sort_interactions <- function(df) {
  df.res <- data.frame()
  for (i in 1:nrow(df)) {
    print(sprintf("#%s/#%s", i, nrow(df)))
    
    chrA <- df[i,1]
    posA <- df[i,2]
    
    chrB <- df[i,3]
    posB <- df[i,4]
    
    ### Ordering
    idx.order <- order(c(chrA, chrB))
    
    ### Chr
    chrA.order <- c(chrA, chrB)[idx.order][1] # take first element of order
    chrB.order <- c(chrA, chrB)[idx.order][2]
    ### Pos
    posA.order <- c(posA, posB)[idx.order][1] # take first element of order
    posB.order <- c(posA, posB)[idx.order][2]
    
    #print(c(chrA.order, posA.order, chrB.order, posB.order))
    #row <- matrix(c(chrA.order, posA.order, chrB.order, posB.order, interactionID))
    
    row <- data.frame(chrA=chrA.order, posA=posA.order, chrB=chrB.order, posB=posB.order)
    df.res <- rbind.fill(df.res, row)
  }
  return(df.res)
}
