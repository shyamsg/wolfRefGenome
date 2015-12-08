##############################################################
# File: estSFS.dosage.R                                      #
# Author: Shyam Gopalakrishnan                               #
# Date: 1st December 2015                                    #
# Description: Script to compute the SFS from the dosage     #
# files generated from the vcf files.                        #
##############################################################

args = commandArgs(trailingOnly=TRUE)
dosageFile = args[1]
memFile = args[2]
outfile = args[3]

maf = function(curDosage) {
    af = mean(curDosage, na.rm=T)/2.0
    return(min(af, 1-af))
}

## Read the dosage file.
dosages = read.table(dosageFile, header=T)

## Read the memfile
clusterMem = read.table(memFile, as.is=TRUE, header=TRUE)

## Compute maf
dog.dosages = dosages[, which(clusterMem$species == "Dog")]
wolf.dosages = dosages[, which(clusterMem$species == "Wolf")]

dog.mafs = apply(dog.dosages, 1, maf)
wolf.mafs = apply(wolf.dosages, 1, maf)

pdf(outfile, width=14, height=7)
layout(matrix(1:2, nr=1))
hist(dog.mafs, breaks=seq(0,0.5,0.05), xlim=c(0,0.5), main="SFS in dogs", xlab="MAF", col="firebrick")
hist(wolf.mafs, breaks=seq(0,0.5,0.05), xlim=c(0,0.5), main="SFS in wolves", xlab="MAF", col="forestgreen")
dev.off()