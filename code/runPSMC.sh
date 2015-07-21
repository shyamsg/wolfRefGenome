#! /bin/bash

# File: runPSMC.sh
# Author: Shyam Gopalakrishnan
# Date: 15th July 2015
# This script is the pipeline to run psmc on 
# all the samples in the data folder

# Generate the psmcfa file for the sample from the single sample all sites vcf data
for i in /home/shyam/projects/wolfRefGenome/data/*/*vcf.gz; do 
    /home/shyam/projects/wolfRefGenome/code/vcf2psmc.sh -m 6 -w 100 -b -v $i >& /home/shyam/projects/wolfRefGenome/logs/`basename $i .vcf.gz`.vcf2psmc.log & 
done

# Run psmc for the samples that have the psmcfa file
for i in /home/shyam/projects/wolfRefGenome/data/*/*psmcfa; do 
    /home/shyam/install/packages/psmc/psmc -N20 -t10 -r5 -p "1*6+58*1" -o /home/shyam/projects/wolfRefGenome/analysis/psmc/$(basename $i fa) $i >& /home/shyam/projects/wolfRefGenome/logs/$(basename $i fa).log &
done
