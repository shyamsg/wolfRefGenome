#! /bin/bash

# File: runPSMC.sh
# Author: Shyam Gopalakrishnan
# Date: 15th July 2015
# This script is the pipeline to run psmc on 
# all the samples in the data folder

module load bcftools/1.2

MINDEPTH=$1
WINDOWSIZE=$2

# Generate the psmcfa file for the sample from the single sample all sites vcf data
for i in /home/shyam/projects/wolfRefGenome/data/dogRef/*.depthdist.vcf.gz; do
    rootname=`basename $i .vcf.gz`
    if [ ! -e /home/shyam/projects/wolfRefGenome/data/dogRef/$rootname.win$WINDOWSIZE.miss80.mind$MINDEPTH.psmcfa ]; then
	/home/shyam/projects/wolfRefGenome/code/vcf2psmc.sh -m $MINDEPTH -w $WINDOWSIZE -b -v $i >& /home/shyam/projects/wolfRefGenome/logs/$rootname.vcf2psmc.log &
    fi
done

for i in /home/shyam/projects/wolfRefGenome/data/wolfRef/*depthdist.vcf.gz; do
    rootname=`basename $i .vcf.gz`
    if [ ! -e /home/shyam/projects/wolfRefGenome/data/wolfRef/$rootname.win$WINDOWSIZE.miss80.mind$MINDEPTH.psmcfa ]; then
	/home/shyam/projects/wolfRefGenome/code/vcf2psmc.sh -m $MINDEPTH -w $WINDOWSIZE -b -v $i >& /home/shyam/projects/wolfRefGenome/logs/$rootname.vcf2psmc.log &
    fi
done

# wait for the psmcfa files to be generated
wait
echo "Generated psmcfa files."

# Run psmc for the samples that have the psmcfa file
for i in /home/shyam/projects/wolfRefGenome/data/dogRef/*psmcfa; do
    rootname=`basename $i fa`
    if [ ! -e /home/shyam/projects/wolfRefGenome/analysis/psmc/$rootname ]; then
	/home/shyam/install/packages/psmc/psmc -N20 -t10 -r5 -p "1*6+58*1" -o /home/shyam/projects/wolfRefGenome/analysis/psmc/$rootname $i >& /home/shyam/projects/wolfRefGenome/logs/$rootname.log &
    fi
done

for i in /home/shyam/projects/wolfRefGenome/data/wolfRef/*psmcfa; do
    rootname=`basename $i fa`
    if [ ! -e /home/shyam/projects/wolfRefGenome/analysis/psmc/$rootname ]; then
	/home/shyam/install/packages/psmc/psmc -N20 -t10 -r5 -p "1*6+58*1" -o /home/shyam/projects/wolfRefGenome/analysis/psmc/$rootname $i >& /home/shyam/projects/wolfRefGenome/logs/$rootname.log &
    fi
done
