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

## Dog ref mapping vcfs from paleomix - done by Sama
## The files are named bgz but Sama assures me that
## they are gz, so renaming the link.
cd $DATA_HOME/dogRef
## Novembre data
for vcf in /disk/lemaitre/data/joseas/Wolf/Mikkel/Novembre_data/VCFsDog/results/canines_ref_dog/genotypes/*filtered.vcf.bgz; do
    if [ ! -e $vcf ]; then
	ln -s $vcf $(basename $vcf | sed 's/bgz/gz/')
    fi
done
## Wang data
for vcf in /home/mischu/data/projects/2014-03_ancient_dogs/phylogeny/results/canines/genotypes/W*ed.vcf.bgz; do
    if [ ! -e $vcf ]; then
	ln -s $vcf $(basename $vcf | sed 's/bgz/gz/')
    fi
done
## Zhang data
for vcf in /home/joseas/data/Wolf/Mikkel/Zhang_data/VCFsDog/results/canines_ref_dog/genotypes/*filtered.vcf.bgz; do
    if [ ! -e $vcf ]; then
	ln -s $vcf $(basename $vcf | sed 's/Chinese/Zhang/' | sed 's/bgz/gz/').
    fi
done

