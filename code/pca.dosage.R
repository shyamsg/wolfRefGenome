##############################################################
# File: pca.dosage.R                                         #
# Author: Shyam Gopalakrishnan                               #
# Date: 1st December 2015                                    #
# Description: Script to perform PCA using the dosage files  #
# generated from the vcf files.                              #
##############################################################

args = commandArgs(trailingOnly=TRUE)
dosageFile = args[1]
memFile = args[2]
minMaf = as.numeric(args[3])
outfile = args[4]

maf = function(curDosage) {
    af = mean(curDosage, na.rm=T)/2.0
    return(min(af, 1-af))
}

isVariant = function(curDosage) {
    return(var(curDosage, na.rm=T) > 0)
}

## Read the dosage file.
dosages = read.table(dosageFile, header=T)

## Filter out based on maf
dosages = na.omit(dosages)
mafs = apply(dosages, 1, maf)
dosages = t(dosages[which(mafs > minMaf),])
which.variant = apply(dosages, 2, isVariant)
dosages = dosages[, which.variant]

## Read the cluster file
clusterMem = read.table(memFile, header=T)
cols = as.character(clusterMem$species)
cols[cols=="Dog"] = "firebrick"
cols[cols=="Wolf"] = "forestgreen"
cols[cols=="Jackal"] = "darkgray"
pchs = as.numeric(clusterMem$dataset)+15

## Perform PCA on the dosages
pca.dosage = prcomp(dosages, center=TRUE, scale=TRUE, retx=TRUE)
pca.vars = pca.dosage$sdev**2
pca.vars = round(pca.vars*100.0/sum(pca.vars), 2)

## Plot the PCs
pdf(file=outfile, width=9, height=12)
layout(matrix(1:12, nr=4, byrow=TRUE))
for (index in 1:11) {
    p1 = 2*index - 1
    p2 = p1 + 1
    plot(pca.dosage$x[,p1], pca.dosage$x[,p2], xlab=paste0("PC",p1, " (",pca.vars[p1], ")"),
         ylab=paste0("PC",p2, " (",pca.vars[p2], ")"),
	 pch=pchs, col=cols)
    legend("topright", legend=c("Dog","Wolf","Jackal"), col=c("firebrick","forestgreen","darkgray"), pch=15)
    legend("bottomright", legend=levels(clusterMem$dataset), pch=16:18, col=2)
}
dev.off()
