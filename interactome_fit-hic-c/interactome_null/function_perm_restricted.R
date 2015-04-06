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
