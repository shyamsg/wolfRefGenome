######################################################
# File: pca.ngscovar.R                               #
# Author: Shyam Gopalakrishnan                       #
# Description: Plot the pca of a given ngscovar made #
# var-covar matrix using the given cluster file.     #
######################################################

args = commandArgs(trailingOnly=TRUE)
covarfile = args[1]
memfile = args[2]
outpdf = args[3]

covar = as.matrix(read.table(covarfile))
pca.covar = eigen(covar)
pca.covar$values = pca.covar$values/sum(pca.covar$values)
pca.covar$values = round(pca.covar$values*100, 2)

## Get colors and shapes from the mems file
mems = read.table(memfile, header=T, as.is=T)

cols = as.numeric(as.factor(mems$species))
pchs = as.numeric(as.factor(mems$dataset))+15

## Plot the PCA plot
pdf(file=outpdf)
plot(pca.covar$vectors[,1], pca.covar$vectors[,2], col=cols, pch=pchs, xlab=paste0("PC 1 (", pca.covar$values[1], "%)"), ylab=paste0("PC 2 (", pca.covar$values[2],"%)"), cex=0.7)
legend("topright", legend=levels(as.factor(mems$species)), pch=16, col=1:3)
legend("bottomright", legend=levels(as.factor(mems$dataset)), pch=16:18, col=1)
plot(pca.covar$vectors[,3], pca.covar$vectors[,4], col=cols, pch=pchs, xlab=paste0("PC 3 (", pca.covar$values[3], "%)"), ylab=paste0("PC 4 (", pca.covar$values[4],"%)"), cex=0.7)
legend("topright", legend=levels(as.factor(mems$species)), pch=16, col=1:3)
legend("bottomright", legend=levels(as.factor(mems$dataset)), pch=16:18, col=1)
dev.off()