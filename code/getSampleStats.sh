#! /bin/bash

# File: getSampleStats.sh
# Author: Shyam Gopalakrishnan
# Date: 21 July 2015
# This script gets the individual level summary statistics
# from all the single sample vcfs. It is simply a wrapper
# for the vcftools options. 

# Generate tabix indexes for the vcf files
for invcf in /home/shyam/projects/wolfRefGenome/data/*/*vcf.gz; do 
    if [ ! -e $invcf.tbi ]; then
        tabix -p vcf $invcf &
    fi
done

VCFTOOLS='/usr/local/src/vcftools/0.1.12b/bin/vcftools'
# Compute the heterozygosity, missingness and depth statistics for each sample
for invcf in /home/shyam/projects/wolfRefGenome/data/*/*vcf.gz; do
    plot_sampleStats.py -v $invcf -o /home/shyam/projects/wolfRefGenome/analysis/$(basename $invcf .vcf.gz).pdf
done

export PERL5LIB=/usr/local/src/vcftools/0.1.12b/perl
VCFMERGE='/usr/local/src/vcftools/0.1.12b/perl/vcf-merge'
# Combined the dog ref and the wolf ref into two multisample vcfs.
for direc in dogRef wolfRef; do 
    $VCFMERGE /home/shyam/projects/wolfRefGenome/data/$direc/*vcf.gz | bgzip -c > /home/shyam/projects/wolfRefGenome/data/$direc/all_$direc.vcf.gz &
done
