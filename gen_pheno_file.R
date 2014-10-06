rm(list=ls())


wd = '/Users/pascaltimshel/Dropbox/EGCUT_DATA/Expression_related_docs_pascal'
setwd(wd)

file.in = 'ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.with_probes246.extract.txt'
df.in <- read.delim(file.in, row.names=1, check.names=F)
df.in.t <- t(df.in)

file.out <- 'ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.with_probes246.extract.transpose.txt'
#write.csv(df.in.t, file=file.out)
write.table(df.in.t, file=file.out, quote=F, sep="\t")

?read.table
