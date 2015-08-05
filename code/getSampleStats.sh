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

# Modify the vcf file to remove the null indels
# meaning that indels where it is called an indel
# but the alternate allele is not different from
# reference and the PL at the end of the line has
# only one value, is removed from the file.
for invcf in /home/shyam/projects/wolfRefGenome/data/*/*vcf.gz; do 
    if [ ! -e ${invcf/.vcf.gz/nullRemoved.vcf.gz} ]; then
	echo /home/shyam/projects/wolfRefGenome/code/removeNullIndels.py -v $invcf -o ${invcf/.vcf.gz/.nullRemoved.vcf.gz} 
    fi 
done


export PERL5LIB=/usr/local/src/vcftools/0.1.12b/perl
VCFTOOLS='/usr/local/src/vcftools/0.1.12b/bin/vcftools'
VCFMERGE='/usr/local/src/vcftools/0.1.12b/perl/vcf-merge'

# For each sample, we can filter the vcf file and create a new
# file that does NOT have indels, snps that are close to indels, 
# snps that are close to each other
for invcf in /home/shyam/projects/wolfRefGenome/data/*/Wang*.vcf.gz; do 
    $VCFTOOLS --gzvcf $invcf --stdout --recode --remove-indels --recode-INFO-all | bgzip -c > ${invcf/.vcf.gz/.noindel.vcf.gz} &
    $VCFTOOLS --gzvcf $invcf --stdout --recode --keep-only-indels --recode-INFO-all | bgzip -c > ${invcf/.vcf.gz/.indels.vcf.gz} &
done


# Combined the dog ref and the wolf ref into two multisample vcfs.
for direc in dogRef; do 
    $VCFMERGE -d /home/shyam/projects/wolfRefGenome/data/$direc/W*vcf.gz | bgzip -c > /home/shyam/projects/wolfRefGenome/data/$direc/all_$direc.vcf.gz &
done
