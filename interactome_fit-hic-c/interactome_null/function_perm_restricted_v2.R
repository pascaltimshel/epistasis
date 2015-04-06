###################### Null - restrictions v.2 #################

perm_restricted <- function(n_perm, df.interaction_table) {
  #################################### DEPENDENCIES ####################################
  library(dplyr)
  
  time.function.start <- proc.time()
  #################################### DESCRIPTION ####################################
  ### Input
  # n_perm:                 number of permutations to return
  # df.interaction_table :  a data frame containing the columns "chrA", "posA", "chrB", "posB". Other columns will be ignored. 
  #                         this data frame XXX x will be used to avoid any "positives".
  
  ### Output
  # r.list: contains two data frames: df.perm, df.perm.stats
  
  #################################### FUNCTION ####################################
  
  ################## Parameters ##################
  param.null.inter_chromosomal <- TRUE      # bool. if TRUE, then the null will only contain inter chromosomal interactions
  param.distance_restriction <- TRUE        # bool. if TRUE, then the null will enforce a distance restriction given by "param.distance_restriction.value", to ensure that the partnerB and "permuted partnerB" is seperated by a minimum distance. Without this option enabled, the null and hi-c could potentially share SNP-pairs
  param.distance_restriction.value <- 1e5  # integer (bp). Set the minimal distance between any interactions 
                                            # *OBS* this parameter is only active if "param.distance_restriction" is TRUE
                                            # 1e5: 100 kb
                                            # 1e6 1 Mb
  
  #param.null.same_chromosome_pair <- FALSE  # *NOT IMPLEMENTED YET*. Determined if the Hi-C and null should have the chromosome pair forming the interaction.
  
  ################## Initialyzation ##################
  ### For debugging - running for-loop manually
  #n_perm <- 10
  
  ### Set seed
  set.seed(1) # Important to set seed for REPRODUCIBLE results
  
  ### number of interactions
  n_interactions <- nrow(df.interaction_table)
  
  ### ADDING NEW COLUMN: Setting a unambiguosly identifer for interactions
  df.interaction_table$interaction_identifer <- with(df.interaction_table, paste0(chrA,":", posA, "_", chrB,":", posB)) # e.g. "1:171570144_9:118763863"
  # ^ this only gives a  because the interaction table have SORTED partnerA and partnerB to *UNAMBIGUOUSLY* represent an interaction
  ### ADDING NEW COLUMN: index
  df.interaction_table$index <- seq(1, n_interactions)
  
  ### Initialyzing df/matrix that will contain null (permuted indecies)
  # rows: LENGTH=n_interactions
  # columns: LENGTH=n_perm. columns contains "null indexes". integer; can take on values 1..n_interactions.
  df.perm <- data.frame(matrix(NA, ncol=n_perm, nrow=n_interactions))
  colnames(df.perm) <- paste0("null_", seq(n_perm))
  rownames(df.perm) <- paste0("interaction_", seq(n_interactions)) # data frame could also inherit rownames directly from "df.interaction_table"

  ### Stat data frame
  df.perm.stats <- data.frame(matrix(NA, ncol=0, nrow=n_interactions)) # empty data frame

  ################## DATA Validation checks ##################
  ### Checking for duplicates - stop execution if duplicates are detected
  stopifnot(!anyDuplicated(subset(df.interaction_table, select = c("chrA", "posA", "chrB", "posB"))))
  
 ################## Main loop - *NEW VERSION* ##################
  for (interaction_no in (1:n_interactions)) {
    ### ^^OBS: it is important loop "over rows". This allows us to make sure that there are no duplicate index per row.
    
    ### Display mesage
    print(sprintf("#%s/#%s", interaction_no, n_interactions))
    
    ######## Defining variables ########
    ### Extracting variables
    now.interaction_identifer <- df.interaction_table[interaction_no, "interaction_identifer"]
    now.partnerA.chr <- df.interaction_table[interaction_no, "chrA"] # needed for inter-chromosomal
    now.partnerA.pos <- df.interaction_table[interaction_no, "posA"]
    now.partnerB.chr <- df.interaction_table[interaction_no, "chrB"] # needed for distance criteria
    now.partnerB.pos <- df.interaction_table[interaction_no, "posB"] # needed for distance criteria
    
    ### "null identifier" - generate interaction identifiers for current partnerA and all partnerBs
    # used for Selection step #1: exclude existing interactions
    null.interaction_identifer <- paste0(now.partnerA.chr, ":", now.partnerA.pos, "_", df.interaction_table$chrB, ":", df.interaction_table$posB) # vector
    
    ### Stats
    df.perm.stats$interaction_identifer[interaction_no] <- now.interaction_identifer
    df.perm.stats$n_pool.initial[interaction_no] <- n_interactions
    
    ######## Selection - step #1: exclude any existing interactions from the interaction table ########
    # Select partnerBs (interactions), that does *NOT* form interaction pairs in the interaction table
    pool_valid_partners <- df.interaction_table %>%
      filter(interaction_identifer != null.interaction_identifer)
    
    now.n_pool.s1 <- nrow(pool_valid_partners)
    
    ### Saving stats
    stat.n_pool.s1.shrinkage <- n_interactions - now.n_pool.s1
    df.perm.stats$s1.shrinkage[interaction_no] <- stat.n_pool.s1.shrinkage
    
#     if (now.n_pool.s1 < n_interactions-1 ) { # inform the user, if our pool has shrunken more than one interaction. NB: we expect the pool to decrease by one, because the interaction itself is excluded.
#       print(sprintf("no #%s/#%s | Selection step #1 | pool_valid_partners = %s", interaction_no, n_interactions, now.n_pool.s1 ))
#     }
    
    ######## Selection - step #2: exclude INTRA-chromosomal interactions ########
    # select partnerB that are *NOT* on same chromosome as partnerA
    if (param.null.inter_chromosomal) {
      pool_valid_partners <- pool_valid_partners %>%
        filter(chrA != now.partnerA.chr)
    }
    now.n_pool.s2 <- nrow(pool_valid_partners)
    
    ### Saving stats
    stat.n_pool.s2.shrinkage <- now.n_pool.s1 - now.n_pool.s2
    df.perm.stats$s2.shrinkage[interaction_no] <- stat.n_pool.s2.shrinkage
    
    ######## Selection - step #3: distance criteria ########
    if (param.distance_restriction) {
      bool.valid_distance.chr <- with(pool_valid_partners, (chrB != now.partnerB.chr))
      bool.valid_distance.upstream <- with(pool_valid_partners, (chrB == now.partnerB.chr) & (posB <= now.partnerB.pos - param.distance_restriction.value)) # upstream: smaller numbers
      bool.valid_distance.downstream <- with(pool_valid_partners, (chrB == now.partnerB.chr) & (posB >= now.partnerB.pos + param.distance_restriction.value)) # downstream: larger numbers
      ## ^^ explanation of "(chrB==now.partnerB.chr) & ...": distance does not matter if now.partnerB.chr is identical to the null
      pool_valid_partners <- pool_valid_partners %>%
        filter(bool.valid_distance.chr | bool.valid_distance.upstream | bool.valid_distance.downstream)
    }
    
    now.n_pool.s3 <- nrow(pool_valid_partners)
    ### Saving stats
    stat.n_pool.s3.shrinkage <- now.n_pool.s2 - now.n_pool.s3
    df.perm.stats$s3.shrinkage[interaction_no] <- stat.n_pool.s3.shrinkage
    #if (now.n_pool.s3 < now.n_pool.s2) { # inform that our pool have decreased in size
      #print(sprintf("no #%s/#%s | Selection step #3 | pool_valid_partners = %s", interaction_no, n_interactions, now.n_pool.s3 ))
    #}
    

    n_pool.final <- now.n_pool.s3
    df.perm.stats$n_pool.final[interaction_no] <- n_pool.final

    ######## Sampling ########
    if (now.n_pool.s3 < n_perm) { # the sampling pool is too small to sample without replacement.
      # TO/*OBS*: the sample size could be enlarged by allowing sampling from partnerAs as well.
      #warning(sprintf("#%s/#%s | %s | Sample pool size (%s) is less than n_perm. Sampling WITH replacement. The null will contain duplicates", interaction_no, n_interactions, now.interaction_identifer, now.n_pool.final))
      null.idx <- sample(pool_valid_partners$index, n_perm, replace=TRUE)
      df.perm.stats$sampling_with_replacement[interaction_no] <- TRUE
    } else {
      # sample the number of null index you need WITHOUT replacement --> no duplicates.
      null.idx <- sample(pool_valid_partners$index, n_perm, replace=FALSE)
      df.perm.stats$sampling_with_replacement[interaction_no] <- FALSE
    }
    
    ######## Assignment ########
    # assign to row in df/matrix
    df.perm[interaction_no, ] <- null.idx
  
  } ### END FOR LOOP
 
  ################## Information ##################
  print(sprintf("Number of sampling with replacement: %s", sum(df.perm.stats$sampling_with_replacement)))
  print(sprintf("Minimum sampling pool size: %s", min(df.perm.stats$n_pool.final)))
  time.function.elapsed <- proc.time()[3] - time.function.start[3]
  print(sprintf("Done in %.2f sec [%.2f min]", time.function.elapsed, time.function.elapsed/60))

  ################## Return values ##################
  r.list <- list(df.perm, df.perm.stats)
  return(r.list)
}





#     ### generating random numbers between 1 and "n_interactions"
#     sample_pool <- sample(n_interactions)
#     ### removing the interaction_no from the available indecies. 
#     # this is because the current interaction_no corresponds to an index for the Hi-C data (observed)
#     sample_pool <- sample_pool[!sample_pool %in% interaction_no]
#     ### sample the number of null index you need WITHOUT replacement --> no duplicates.
#     null_idx <- sample(sample_pool, size=n_perm)
#     ### assign to row in matrix
#     matrix_perm_idx[interaction_no, ] <- null_idx

################## Main loop - OLD VERSION *INCOMPLETE* ##################
#   for (i in 1:n_perm) {
#     n_perm_attempts <- 1
#     flag_make_new_permutation <- TRUE
#     while (flag_make_new_permutation) {
#       #perm <- integer(length = n_interactions) # initialyze with 0 values
#       perm <- vector(length = n_interactions) # initialyze with FALSE values
#       for (interaction_no in (1:n_interactions)) {
#         perm[interaction_no] <- XXXX
#       #perm <- sample.int(n_interactions) # this will sample integers in the range 1..."n_interactions".
#       perm_val <- x[perm] # now reorder x according to the permuted index in perm.
#       
#       bool_identical_positions <- x == perm_val
#       print(sprintf("i=%s | Attempt=%s | %s identical positions in x and perm", i, n_perm_attempts, sum(bool_identical_positions)))
#       flag_make_new_permutation <- any(bool_identical_positions) # value will be TRUE if any positions are identical    
#       n_perm_attempts <- n_perm_attempts + 1
#     }
#     print(sprintf("i=%s | Ok, found perm", i))
#     df.perm[,i] <- perm
#     df.perm.stats[i,] <- n_perm_attempts
#   }
#   
#   df.perm.stats$error <- sapply(seq_along(df.perm), function(i) {sum(df.perm[,i]==x)/dim(df.perm)[1]})
#   
