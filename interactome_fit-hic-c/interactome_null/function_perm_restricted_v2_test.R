n_interactions <- length(x)

### Initialyzing matrix with index
# rows: n_interactions
# columns: n_perm. columns contains "null indexes". integer; can take on values 1..n_interactions.
matrix_perm_idx <- matrix(NA, ncol=n_perm, nrow=n_interactions)
for (interaction_no in (1:n_interactions)) {
  ### ^^OBS: it is important loop "over rows". This allows us to make sure that there are no duplicate index per row.
  
  ### generating random numbers between 1 and "n_interactions"
  sample_pool <- sample(n_interactions)
  ### removing the interaction_no from the available indecies. 
  # this is because the current interaction_no corresponds to an index for the Hi-C data (observed)
  sample_pool <- sample_pool[!sample_pool %in% interaction_no]
  ### sample the number of null index you need WITHOUT replacement --> no duplicates.
  null_idx <- sample(sample_pool, size=n_perm)
  ### assign to row in matrix
  matrix_perm_idx[interaction_no, ] <- null_idx
}

matrix_perm_idx

###############################################################
################## Initializing list_available_idx - *NOT USED* ##################
list_available_idx <- list()
for (interaction_no in (1:n_interactions)) {
  available_idx <- sample(n_interactions) # generating random numbers between 1 and "n_interactions"
  ### removing the interaction_no from the available indecies. 
  # this is because the current interaction_no corresponds to an index for the Hi-C data (observed)
  available_idx <- available_idx[!available_idx %in% interaction_no]
  list_available_idx[[interaction_no]] <- available_idx
}
list_available_idx[[1]] ==1 
