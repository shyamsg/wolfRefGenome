Wolf Reference Genome Project
=============================

This document describes the processing done for the reference
genome project. It describes the variant calling processing 
and the subsequent analyses that need to be done to include 
in the de novo reference genome project. 

Directory structure
-------------------
* code

This directory contains all the code developed for this project. 

* data

The data is stored in the this directory in two sub-directories.
One for the reads mapped to the wolf reference genome and another
for the reads mapped to the dof reference genome.

* analysis

Each kind of analysis is stored in a separate sub-directory in this
directory.

Data processing
---------------
Each sample is processed using the same processing pipeline. 
The samples are processed, each independently, to get the 
variants in each sample. The variants from each sample are 
combined to get the list of variants in all the samples. 

Filtering
---------
The variants are filtered to ensure that the variants retained 
are the ones that are reliable. The filters might involve the 
identification in multiple samples, the depth of coverage in 
each sample and other such statistics.


Analyses
--------
* PCA
For both the dog and wolf reference genome, we use only the filtered
markers to perform PCA. The PCA should be fairly straightforward.
* Inbreeding estimation
Since we have only 1-3 samples per population, it might not be the best 
idea to estimate the inbreeding in each sample using the allele frequency
estimates. Also, since the allele frequency estimates itself are going to 
be biased due to any inbreeding, it might make even less sense. It might 
be better to use only the sample heterozygosity using sites that we are 
confident of in each samples separately.
* PSMC
This is closely connected to the sample heterozygosity. We can process the 
data to get the information that we need. We need to develop a method to 
convert the all sites vcf to the input format for PSMC. The filtering 
needs to be thought over a little bit. Need to remove (or set to N) sites
that are not super confidently called. Also, how to translate bad sites to 
bad windows?
* Phylogeny and ABBA-BABA tests
The phylogeny can be constructed using a neighbor joining tree. Then ABBA-
BABA tests can be conducted using the golden jackal as the outgroup. It should
be noted that the violation of D != 0 in the ABBA-BABA can mean things in the 
case where the outgroup is not necessarily a proper outgroup. 
