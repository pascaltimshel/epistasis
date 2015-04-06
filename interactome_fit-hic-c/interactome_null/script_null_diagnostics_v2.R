############### SYNOPSIS ###################
# class: function
# description: plotting diagnositcs for the null generation

# dependencies: df.perm and df.perm.stats NEEDS to be defined

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
file.plot <- file.path(path.out.diagnostics, basename.plot)
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
file.plot <- file.path(path.out.diagnostics, basename.plot)
file.plot
ggsave( file=file.plot )

############ Histogram: number of duplicated interactions across all 'experiments'
df.perm.dup <- data.frame(n_duplicates=apply(df.perm, 1, function(row) {sum(duplicated(row))})) # finding the number of indentical null for each interaction
# ggplot
ggplot(df.perm.dup, aes(x=n_duplicates)) + geom_histogram() + labs(title="Distribution of number of duplicated interactions", x="Number of duplicated interactions", y="Frequency")
# save
basename.plot <- sprintf("diagnostics_null_number_of_duplicated_interactions_%s_nperm_%s_q_%s.pdf", hic_cell_type, n_perm, p.q.threshold)
file.plot <- file.path(path.out.diagnostics, basename.plot)
file.plot
ggsave( file=file.plot )

############ Histogram: colSums [DOES NOT MATTER]
colsum.df.perm <- colSums(df.perm)
qplot(colsum.df.perm)