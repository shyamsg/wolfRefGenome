#! /bin/bash 

# File: vcf2psmc.sh
# Author: Shyam Gopalakrishnan
# Date: 11th July 2015
# Script to convert a single sample all sites vcf to 
# a psmcfa input file. 

if [ $HOSTNAME == 'kraken' ]; then
    VCFUTILS='/usr/local/src/samtools/0.1.18/bcftools/vcfutils.pl'
    VCFTOOLS
elif [ $HOSTNAME == 'lemaitre' ]; then
    VCFUTILS='vcfutils.pl'
else 
    VCFUTILS='vcfutils.pl'
fi

infile=""
infilenums=0
minDepth=5
window=100
force_bgzip=0

if [ $# -lt 1 ]; then 
    echo "Usage: vcf2psmc.sh [-m|--minDepth minimumDepth] [-w|--window windowSize] [-b|--bgzip] -v in.vcf|in.vcf.gz|in.vcf.bgz"
    exit 1
fi
while [[ $# > 0 ]]; do
    key="$1"
    case $key in
	-v|--vcf)
	    infile="$2"
	    let infilenums=infilenums+1
	    shift # past argument
	    ;;
	-m|--minDepth)
	    minDepth="$2"
	    shift # past argument
	    ;;
	-w|--window)
	    window="$2"
	    shift # past argument
	    ;;
	-b|--bgzip)
	    force_bgzip=1
	    ;;
	-h|--help)
	    echo "Usage: vcf2psmc.sh [-m|--minDepth minimumDepth] [-w|--window windowSize] [-b|--bgzip] -v in.vcf|in.vcf.gz|in.vcf.bgz"
	    exit 0
	    ;;
	*)
            # unknown option
	    echo "Unknown option: $key"
	    echo "Usage: vcf2psmc.sh [-m|--minDepth minimumDepth] [-w|--window windowSize] [-b|--bgzip] -v in.vcf|in.vcf.gz|in.vcf.bgz"
	    exit 1
	    ;;
    esac
    shift # past argument or value
done

if [ "x$infile" == "x" ]; then 
    echo "A single input file is required. Please use the -v option to input one."
    exit 1
elif [ $infilenums -ne 1 ]; then
    echo "The -v option should be specified only once."
    exit 1
fi

echo "Starting the processing for $infile."
# Check if the file is gz, bgz or just plain old vcf
# Depending on that, process the file to generate the 
# fastq file, which will subsequently processed to get the
# psmcfa file
direc=`dirname $infile`
if [ "x${infile:(-8)}" == "x.vcf.bgz" ]; then
    rootname=${infile:0:(-8)}
    echo "Processing the input bgzipped vcf file."
    bgzip -d -c $infile | $VCFUTILS vcf2fq -d $minDepth | gzip > $rootname.fq.gz
elif [ "x${infile:(-7)}" == "x.vcf.gz" ]; then
    rootname=${infile:0:(-7)}
    if [ $force_bgzip -eq 1 ]; then 
	bgzip -d -c $infile | $VCFUTILS vcf2fq -d $minDepth | gzip > $rootname.fq.gz
    else
	gunzip -c $infile | $VCFUTILS vcf2fq -d $minDepth | gzip > $rootname.fq.gz
    fi
elif [ "x${infile:(-4)}" == "x.vcf" ]; then
    rootname=${infile:0:(-4)}
    $VCFUTILS vcf2fq -d $minDepth $infile | gzip > $rootname.fq.gz
else
    echo "The input file does not have a valid extension. Please correct error."
    exit 2
fi
echo "Finished generating the fastq file from vcf."

# Generate the psmcfa file using the script that 
# accompanies the psmc program
/home/shyam/install/packages/psmc/utils/fq2psmcfa -m 0.80 -s $window $rootname.fq.gz > $rootname.win$window.miss80.mind$minDepth.psmcfa
echo "Finished generating the psmcfa file."
#rm -f $rootname.fq.gz
echo "Cleaned up."
