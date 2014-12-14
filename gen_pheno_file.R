rm(list=ls())


wd = '/Users/pascaltimshel/Dropbox/EGCUT_DATA/Expression_related_docs_pascal'
setwd(wd)

######################## 246 probes #######################
file.in = 'ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.with_probes246.extract.txt'
df.in <- read.delim(file.in, row.names=1, check.names=F)
df.in.t <- t(df.in)

file.out <- 'ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.with_probes246.extract.transpose.txt'
write.table(df.in.t, file=file.out, quote=F, sep="\t")

######################## All probes - adding FID and IID #######################
file.in = 'ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.extract.txt'
df.in <- read.delim(file.in, row.names=1, check.names=F)
df.in.t <- t(df.in)
df.in.t.fid.iid <- data.frame(FID=rownames(df.in.t), IID=rownames(df.in.t), df.in.t, check.names=F)

file.out <- 'ExpressionDataCorrected4GWASPCs.ExpressionData.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.with_GTypes832.extract.transpose.txt'
write.table(df.in.t.fid.iid, file=file.out, quote=F, sep="\t", row.names=F)




