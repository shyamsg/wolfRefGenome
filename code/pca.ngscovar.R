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
pdf(file=outpdf, width=12, height=16)
layout(matrix(1:12, nr=4, byrow=T))
for (pcindex in 1:10) {
    p1 = 2*pcindex-1
    p2 = p1+1
    plot(pca.covar$vectors[,p1], pca.covar$vectors[,p2], col=cols, pch=pchs, xlab=paste0("PC ", p1, " (", pca.covar$values[p1], "%)"), ylab=paste0("PC ", p2, " (", pca.covar$values[p2],"%)"), cex=0.7)
    legend("topright", legend=levels(as.factor(mems$species)), pch=16, col=1:3)
    legend("bottomright", legend=levels(as.factor(mems$dataset)), pch=16:18, col=1)
}
dev.off()