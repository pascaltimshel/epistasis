############### SYNOPSIS ###################
# class: function
# description: UTILITY FUNCTIONS for the null generation.
### Steps
# 1) write df.perm to file

write_df.perm <- function(path.out) {
  ### OUTPUT
  # EXPORT df.perm.concat:  a tab seperated file with indexes for the "partnerB" for the interaction.
  #                         DOWNSTREAM description: interaction pairs for EXPERIMENTS are formed by joining the "partnerA" (index "a") column with the index in the null table (index "b") for the given experiment. Note that for the experiment hic_1 (positives), we always have that index a==b; whereas for null experiments we have a != b
  #                             The original interactions that the NULL INDEX POINT TO (link/refer to) are found in the file, e.g.: "interation_table.fit-hi-c.nosex.interchromosomal.hESC.q_1e-13.txt"
  #                         columns = experiments (hic_1, null_1, null_2,...)
  #                         rows = interaction_no
  #                         *REMARKS #1*: notice that the first column of df.perm.concat is the hic_1 experiment. The indexes in this column are 1, 2, 3, .., n_interactions
  #                         *REMARKS #2*: remember that all experiments for all interaction pairs share the same "partnerA". Only "partnerB" is different.
  

  ### Concatenate data frame 
  df.perm.concat <- data.frame(hic_1=seq(nrow(df.perm)), df.perm)
  
  #####################################
  ### Setting pathname
  file.interaction_table.basename.sans_ext <- file_path_sans_ext(basename(file.interaction_table)) # removing trailing ".txt"
  file.interaction_table.basename.sans_ext
  file.null_table.basename.sans_ext <- sub('interation_table', 'null_table', file.interaction_table.basename.sans_ext) # regex substitution
  #file.null_table.sans_ext <- file.path(dirname(file.interaction_table), file.null_table.basename.sans_ext)
  file.null_table.sans_ext <- file.path(path.out, file.null_table.basename.sans_ext)
  file.null_table <- sprintf("%s.nperm_%s.txt", file.null_table.sans_ext, n_perm)
  file.null_table
  
  ### Write table
  write.table(df.perm.concat, file=file.null_table, col.names=T, row.names=F, quote=F, sep="\t")
}