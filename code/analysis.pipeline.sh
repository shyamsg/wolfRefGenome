#! /bin/bash

##################################################
# File: analysis.pipeline.sh                     #
# Author: Shyam Gopalakrishnan                   #
# Date: 25th September 2015                      #
# Description: This shell script contains all    #
# the code that was run by Shyam, starting from  #
# making the directories to store data, to the   #
# final analysis. Of course, this script was not #
# run all at once, since many of the called      #
# functions and scripts were developed along the #
# way, but this script can _technically_ be used #
# to recreate the entire analysis.               #
##################################################

## Top level variables such as root project directory
## and tool names and locations and loading of any
## modules that might be requires
## The root directory should already exists at this point.
## OR things will not move forward. Of course,
## if someone else is recreating the analysis,
## this needs to be changed :) (STOP trying to write
## into my directory)
PROJECT_HOME=/home/shyam/projects/wolfRefGenome
DATA_HOME=/home/shyam/data/wolfRefGenome

##################################################
#                Directory setup                 #
##################################################
## If data directory and data link do not exist in
## the wolfRefGenome directory then make a link.
if [ ! -e $DATA_HOME ]; then
    mkdir $DATA_HOME
fi
if [ ! -e $PROJECT_HOME/data ]; then
    ln -s $DATA_HOME $PROJECT_HOME/data
fi

## Create mid level data directories and link the
## VCFs created by sama into these directories,
## simultaneously renaming the Chinese* wolves to
## Zhang* wolves in order to keep the author name
## convention used in the other two data sets.
if [ ! -d $DATA_HOME/dogRef ]; then
    mkdir $DATA_HOME/dogRef
fi

if [ ! -d $DATA_HOME/wolfRef ]; then
    mkdir $DATA_HOME/wolfRef
fi

echo "Directory setup complete."

#################################################
#  Data linking from Sama's initial processing  #
#################################################
## Dog ref mapping vcfs from paleomix - done by Sama
## Linking the bgzip of the vcf and the tabix of it
cd $DATA_HOME/dogRef
## Novembre data
for vcf in /disk/lemaitre/data/joseas/Wolf/Mikkel/Novembre_data/VCFsDog/results/canines_ref_dog/genotypes/*filtered.vcf.bgz; do
    if [ ! -e $(basename $vcf) ]; then
	echo "Linking Novembre dog mapped data: $vcf."
	ln -s $vcf .
	ln -s $vcf.tbi .
    fi
done
## Wang data
for vcf in /home/mischu/data/projects/2014-03_ancient_dogs/phylogeny/results/canines/genotypes/W*ed.vcf.bgz; do
    if [ ! -e $(basename $vcf) ]; then
	echo "Linking Wang dog mapped data: $vcf."
	ln -s $vcf .
	ln -s $vcf.tbi .
    fi
done
## Zhang data
for vcf in /home/joseas/data/Wolf/Mikkel/Zhang_data/VCFsDog/results/canines_ref_dog/genotypes/*filtered.vcf.bgz; do
    if [ ! -e $(basename $vcf | sed 's/Chinese/Zhang/') ]; then
	echo "Linking Zhang dog mapped data: $vcf."
	ln -s $vcf $(basename $vcf | sed 's/Chinese/Zhang/')
	ln -s $vcf.tbi $(basename $vcf.tbi | sed 's/Chinese/Zhang/')
    fi
done

## Wolf ref mapping vcfs from paleomix - done by Sama
## Linking the bgzip of the vcf and the tabix of it
cd $DATA_HOME/wolfRef
## Novembre data
for vcf in /disk/lemaitre/data/joseas/Wolf/Mikkel/Novembre_data/VCFsWolf/results/canines_ref_wolf/genotypes/*filtered.vcf.bgz; do
    if [ ! -e $(basename $vcf) ]; then
	echo "Linking Novembre wolf mapped data: $vcf."
	ln -s $vcf .
	ln -s $vcf.tbi .
    fi
done
## Wang data
for vcf in /disk/lemaitre/data/joseas/Wolf/VCFsWolf/results/Wang_canines_ref_wolf/genotypes/*filtered.vcf.bgz; do
    if [ ! -e $(basename $vcf) ]; then
	echo "Linking Wang wolf mapped data: $vcf."
	ln -s $vcf .
	ln -s $vcf.tbi .
    fi
done
## Zhang data
for vcf in /home/joseas/data/Wolf/Mikkel/Zhang_data/VCFsWolf/results/canines_ref_wolf/genotypes/*filtered.vcf.bgz; do
    if [ ! -e $(basename $vcf | sed 's/Chinese/Zhang/') ]; then
	echo "Linking Zhang wolf mapped data: $vcf."
	ln -s $vcf $(basename $vcf | sed 's/Chinese/Zhang/')
	ln -s $vcf.tbi $(basename $vcf.tbi | sed 's/Chinese/Zhang/')
    fi
done

echo "Data linking complete."

##################################################
#      Compute the stats for each vcf file       #
##################################################
## Compute the average depth of each sample using vcftools
## all dog mapped vcfs
module load bcftools/1.2
cd $DATA_HOME/dogRef
for vcf in *.vcf.bgz; do
    if [ ! -e $(basename $vcf .vcf.bgz).stats ]; then
	echo "Computing stats for $vcf."
	bcftools stats -d 0,500,1 $vcf >& $(basename $vcf .vcf.bgz).stats &
    fi
done

## all wolf mapped vcfs
cd $DATA_HOME/wolfRef
for vcf in *.vcf.bgz; do
    if [ ! -e $(basename $vcf .vcf.bgz).stats ]; then
	echo "Computing stats for $vcf."
	bcftools stats -d 0,500,1 $vcf >& $(basename $vcf .vcf.bgz).stats &
    fi
done

## This command makes the shell wait till the processes are done. Only
## here so that the next steps wait
wait

echo "VCF statistics generation complete."

##################################################
#       Filtering of vcf based on depth and      #
#             distance to snps/indels            #
##################################################
## Generate a file with the depth of each sample in the directory
cd $DATA_HOME/dogRef
if [ ! -e dogRef.allSamples.avgdepths ]; then
    for statfile in *.stats; do
	$PROJECT_HOME/code/computeAvgDepthFromStats.py -s $statfile >> dogRef.allSamples.avgdepths
    done
fi

cd $DATA_HOME/wolfRef
if [ ! -e wolfRef.allSamples.avgdepths ]; then
    for statfile in *.stats; do
	$PROJECT_HOME/code/computeAvgDepthFromStats.py -s $statfile >> wolfRef.allSamples.avgdepths
    done
fi

echo "Stats completed."

## Run the script to filter the vcf file to remove
## null indels and filter based on min/max depths.
## Also filter on proximity to indels and snps.
MINQUAL=20
MINDIST=5
MINDEPTH=5
NEIGHBORQUAL=10

cd $DATA_HOME/dogRef
for vcf in *.vcf.bgz; do  
    if [ ! -e $(basename $vcf .vcf.bgz).depthdist.vcf.gz ]; then
	AVGDEPTH=`grep $(basename $vcf .vcf.bgz) dogRef.allSamples.avgdepths | cut -f2 -d " "`
	MAXDEPTH=`bc <<< "scale=0;$AVGDEPTH*2"`
	(python $PROJECT_HOME/code/filterSampleLevel.py -i $vcf -q $MINQUAL -n $NEIGHBORQUAL -d $MINDEPTH -D $MAXDEPTH -v $MINDIST -o $(basename $vcf .vcf.bgz).depthdist.vcf && bgzip $(basename $vcf .vcf.bgz).depthdist.vcf) >& $PROJECT_HOME/logs/$(basename $vcf .vcf.bgz).depthdist.log &
    fi
done

## This command makes sure the shell waits for the filtering processes to be done.
## Next steps should be run only if these processes are done.
wait

cd $DATA_HOME/wolfRef
for vcf in *.vcf.bgz; do  
    if [ ! -e $(basename $vcf .vcf.bgz).depthdist.vcf.gz ]; then
	AVGDEPTH=`grep $(basename $vcf .vcf.bgz) wolfRef.allSamples.avgdepths | cut -f2 -d " "`
	MAXDEPTH=`bc <<< "$AVGDEPTH*2"`
	(python $PROJECT_HOME/code/filterSampleLevel.py -i $vcf -q $MINQUAL -n $NEIGHBORQUAL -d $MINDEPTH -D $MAXDEPTH -v $MINDIST -o $(basename $vcf .vcf.bgz).depthdist.vcf && bgzip $(basename $vcf .vcf.bgz).depthdist.vcf) >& $PROJECT_HOME/logs/$(basename $vcf .vcf.bgz).depthdist.log &
    fi
done

## This command makes sure the shell waits for the filtering processes to be done.
## Next steps should be run only if these processes are done.
wait

## Tabix index generation
cd $DATA_HOME/dogRef
for vcf in *.depthdist.vcf.gz; do  
    if [ ! -e $vcf.tbi ]; then
	tabix -p vcf $vcf &
    fi
done

## Wait for the tabix commands to finish
wait

cd $DATA_HOME/wolfRef
for vcf in *.depthdist.vcf.gz; do  
    if [ ! -e $vcf.tbi ]; then
	tabix -p vcf $vcf &
    fi
done

## Wait for the tabix commands to finish
wait

echo "Per sample vcf filtering complete."

#############################################
#       Merge the vcfs and get list of      #
#              triallelic sites             #
#############################################
## Sadly bcftools does not work dues to malformed vcf files,
## esp. the absence of certain filter tags in the header.

# module load vcftools/0.1.14
# export PERL5LIB=/usr/local/src/vcftools/0.1.14/perl

# cd $DATA_HOME/dogRef
# if [ ! -e allAlignedToDog.snps_indels.vcf.gz ]; then
#     perl /usr/local/src/vcftools/0.1.14/perl/vcf-merge *.depthdist.vcf.gz | bcftools > allAlignedToDog.snps_indels.vcf.gz
    
# fi

# cd $DATA_HOME/wolfRef
# if [ ! -e allAlignedToWolf.snps_indels.vcf.gz ]; then
#     perl /usr/local/src/vcftools/0.1.14/perl/vcf-merge *.depthdist.vcf.gz | vcftools --vcf - --maf 0.01 --remove-indels --stdout | bgzip -c > allAlignedToWolf.snps_indels.vcf.gz 
# fi

# ## wait for bcftools to finish
# wait

# echo "Made the combined snps only file."

#############################################
#     Process the vcf files to make the     #
#      seq files for psmc and run psmc      #
#############################################
#$PROJECT_HOME/code/runPSMC.sh 

## wait for the PSMC processed to finish.
#wait
