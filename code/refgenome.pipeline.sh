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
export PROJECT_HOME=/home/projects/pr_99009/people/shygop/wolves/wolfRefGenome
export DATA_HOME=$PROJECT_HOME/vcfs

##################################################
#                Directory setup                 #
##################################################
## If data directory and data link do not exist in
## the wolfRefGenome directory then make a link.
if [ ! -e $DATA_HOME ]; then
    mkdir $DATA_HOME
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
# cd $DATA_HOME/dogRef
# ## Novembre data
# for vcf in /disk/lemaitre/data/joseas/Wolf/Mikkel/Novembre_data/VCFsDog/results/canines_ref_dog/genotypes/*filtered.vcf.bgz; do
#     if [ ! -e $(basename $vcf) ]; then
# 	echo "Linking Novembre dog mapped data: $vcf."
# 	ln -s $vcf .
# 	ln -s $vcf.tbi .
#     fi
# done
# ## Wang data
# for vcf in /home/mischu/data/projects/2014-03_ancient_dogs/phylogeny/results/canines/genotypes/W*ed.vcf.bgz; do
#     if [ ! -e $(basename $vcf) ]; then
# 	echo "Linking Wang dog mapped data: $vcf."
# 	ln -s $vcf .
# 	ln -s $vcf.tbi .
#     fi
# done
# ## Zhang data
# for vcf in /home/joseas/data/Wolf/Mikkel/Zhang_data/VCFsDog/results/canines_ref_dog/genotypes/*filtered.vcf.bgz; do
#     if [ ! -e $(basename $vcf | sed 's/Chinese/Zhang/') ]; then
# 	echo "Linking Zhang dog mapped data: $vcf."
# 	ln -s $vcf $(basename $vcf | sed 's/Chinese/Zhang/')
# 	ln -s $vcf.tbi $(basename $vcf.tbi | sed 's/Chinese/Zhang/')
#     fi
# done

# ## Wolf ref mapping vcfs from paleomix - done by Sama
# ## Linking the bgzip of the vcf and the tabix of it
# cd $DATA_HOME/wolfRef
# ## Novembre data
# for vcf in /disk/lemaitre/data/joseas/Wolf/Mikkel/Novembre_data/VCFsWolf/results/canines_ref_wolf/genotypes/*filtered.vcf.bgz; do
#     if [ ! -e $(basename $vcf) ]; then
# 	echo "Linking Novembre wolf mapped data: $vcf."
# 	ln -s $vcf .
# 	ln -s $vcf.tbi .
#     fi
# done
# ## Wang data
# for vcf in /disk/lemaitre/data/joseas/Wolf/VCFsWolf/results/Wang_canines_ref_wolf/genotypes/*filtered.vcf.bgz; do
#     if [ ! -e $(basename $vcf) ]; then
# 	echo "Linking Wang wolf mapped data: $vcf."
# 	ln -s $vcf .
# 	ln -s $vcf.tbi .
#     fi
# done
# ## Zhang data
# for vcf in /home/joseas/data/Wolf/Mikkel/Zhang_data/VCFsWolf/results/canines_ref_wolf/genotypes/*filtered.vcf.bgz; do
#     if [ ! -e $(basename $vcf | sed 's/Chinese/Zhang/') ]; then
# 	echo "Linking Zhang wolf mapped data: $vcf."
# 	ln -s $vcf $(basename $vcf | sed 's/Chinese/Zhang/')
# 	ln -s $vcf.tbi $(basename $vcf.tbi | sed 's/Chinese/Zhang/')
#     fi
# done

echo "Data linking complete."

##################################################
#      Compute the stats for each vcf file       #
##################################################
## Compute the average depth of each sample using vcftools
## all dog mapped vcfs
module load bcftools/1.3.1
cd $DATA_HOME/dogRef
echo "Computing stats."
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
# wait

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
## YAY!!. tweaked the earlier filterSampleLevel.py script
## to remove the filter tags, so that we can use bcftools

module load bcftools/1.2

cd $DATA_HOME/dogRef
if [ ! -e allAlignedToDog.triSNP.bed ]; then
    bcftools merge -m both -Ou *depthdist.vcf.gz | bcftools filter -Ov -i 'MAF[0]>0.01 & TYPE="snp" & ((REF!="N" & N_ALT>1) | N_ALT>2 )' | awk '{ OFS="\t"; if(!/^#/){ print $1,$2-1,$2 }}' > allAlignedToDog.triSNP.bed &
fi

cd $DATA_HOME/wolfRef
if [ ! -e allAlignedToWolf.triSNP.bed ]; then
    bcftools merge -m both -Ou *depthdist.vcf.gz | bcftools filter -Ov -i 'MAF[0]>0.01 & TYPE="snp" & ((REF!="N" & N_ALT>1) | N_ALT>2 )' | awk '{ OFS="\t"; if(!/^#/){ print $1,$2-1,$2 }}' > allAlignedToWolf.triSNP.bed &
fi

## wait for bcftools to finish
wait

echo "Made the trialleleic snps bed file."

############################################
#    Filter out the triallelic snp sites   #
#       single individual vcf files.       #
############################################
## Make another set of per sample vcf files
## that have the trialleleic sites removed.

module load bedtools/2.25.0

cd $DATA_HOME/dogRef
for vcf in *.depthdist.vcf.gz; do
    if [ ! -e $(basename $vcf .vcf.gz).notri.vcf.gz ]; then
	zcat ${vcf/.notri/} | head -20000 | grep "^#" > $(basename $vcf .vcf.gz).notri.vcf && \
	    bedtools subtract -a $vcf -b allAlignedToDog.triSNP.bed >> $(basename $vcf .vcf.gz).notri.vcf && \
	    bgzip $(basename $vcf .vcf.gz).notri.vcf &
    fi
done
wait

cd $DATA_HOME/wolfRef
for vcf in *.depthdist.vcf.gz; do
    if [ ! -e $(basename $vcf .vcf.gz).notri.vcf.gz ]; then
	zcat ${vcf/.notri/} | head -20000 | grep "^#" > $(basename $vcf .vcf.gz).notri.vcf && \
	    bedtools subtract -a $vcf -b allAlignedToWolf.triSNP.bed >> $(basename $vcf .vcf.gz).notri.vcf && \
	    bgzip $(basename $vcf .vcf.gz).notri.vcf &
    fi
done
wait

echo "Made tri removed per sample vcf files."

## Tabix index generation
cd $DATA_HOME/dogRef
for vcf in *.depthdist.notri.vcf.gz; do  
    if [ ! -e $vcf.tbi ]; then
	tabix -p vcf $vcf &
    fi
done

cd $DATA_HOME/wolfRef
for vcf in *.depthdist.notri.vcf.gz; do  
    if [ ! -e $vcf.tbi ]; then
	tabix -p vcf $vcf &
    fi
done
wait

echo "Made tabix files for tri-removed vcfs."

#############################################
#     Process the vcf files to make the     #
#      seq files for psmc and run psmc      #
#############################################
## Remember that the runPSMC.sh script needs
## to have the correct filenames. MAKE SURE
## that these are correct before running this
## section of commands.

MINDEPTH=5
WINDOWSIZE=100

$PROJECT_HOME/code/runPSMC.sh $MINDEPTH $WINDOWSIZE 

## wait for the PSMC processed to finish.
wait

echo "Ran psmc."

##############################################
# Generate the vcf file for input pca or mds #
#   Perform pca and mds, first using plink   #
#        and the using custom R code         #
##############################################

## Use bcftools to merge and filter the data to
## the vcf file satisfying a threshold for mds
## or pca. In this case, we will use no tri sites
## and no sites with N as the reference.

cd $DATA_HOME/dogRef
if [ ! -e allAlignedToDog.snps.bcf ]; then
    bcftools merge -m both -Ou *depthdist.notri.vcf.gz | bcftools filter -Ob -i 'MAF[0]>0.01 & TYPE="snp" & REF!="N" & N_ALT==1' > allAlignedToDog.snps.bcf & 
fi

cd $DATA_HOME/wolfRef
if [ ! -e allAlignedToWolf.snps.bcf ]; then
    bcftools merge -m both -Ou *depthdist.notri.vcf.gz | bcftools filter -Ob -i 'MAF[0]>0.01 & TYPE="snp" & REF!="N" & N_ALT==1' > allAlignedToWolf.snps.bcf & 
fi

## Wait for the bcftools things to finish running.
wait
echo "Made bcf files with all snps."

## Filter this file to get final list of sites
## with specified missingness and maf cutoffs.

module load vcftools/0.1.14
module load bcftools/1.2

MISSING=0.75 ## maximum per site missingness is 25%
MAF=0.05 ## minimum maf
MINQ=30 ## minimum quality
MINGQ=20 ## minimum genotype quality

cd $DATA_HOME/dogRef
if [ ! -e allAlignedToDog.miss$MISSING.maf$MAF.minq$MINQ.mingq$MINGQ.forPCA.bcf ]; then
    bcftools view allAlignedToDog.snps.bcf | vcftools --vcf - --stdout --recode-INFO-all --recode --maf $MAF --max-missing $MISSING --minQ $MINQ --minGQ $MINGQ | sed 's/-1\/-1/.\/./g' | bcftools view -Ob > allAlignedToDog.miss$MISSING.maf$MAF.minq$MINQ.mingq$MINGQ.forPCA.bcf &
fi

cd $DATA_HOME/wolfRef
if [ ! -e allAlignedToWolf.miss$MISSING.maf$MAF.minq$MINQ.mingq$MINGQ.forPCA.bcf ]; then
    bcftools view allAlignedToWolf.snps.bcf | vcftools --vcf - --stdout --recode-INFO-all --recode --maf $MAF --max-missing $MISSING --minQ $MINQ --minGQ $MINGQ | sed 's/-1\/-1/.\/./g' | bcftools view -Ob > allAlignedToWolf.miss$MISSING.maf$MAF.minq$MINQ.mingq$MINGQ.forPCA.bcf &
fi

wait
echo "Made bcf files for pca."

###########################################
# PCA/MDS using the bcf generated in the  #
# step above. First use ngsCovar to get   #
# the covar matrix. Then use R to do pca. #
# Cannot use plink because it cannot take #
# the scaffold genome, since it contains  #
# non standard names for chromosomes.     #
# We also make the dosage file for each of#
# these vcf files and use them for pca.   #
###########################################

DOGFAI='/home/joseas/data/Wolf/RefGenome/canFam31_nucl.fasta.fai'
WOLFFAI='/home/joseas/data/Wolf/VCFsWolf/data/prefixes/L.Dalen_14_wolf.scf.fasta.fai'
NGSCOVAR=/home/fgvieira/data/appz/ngsTools/ngsPopGen/ngsCovar
BCFROOT=$(dirname `which bcftools`)
PLINK19="plink1.9"
VCF2DOSAGE="/home/shyam/data/projects/SantasHelpers/pythonScripts/vcf2Dosage.py"
export LD_LIBRARY_PATH=`dirname $BCFROOT`/htslib-1.2.1/
export BCFTOOLS_PLUGINS=`dirname $BCFROOT`/plugins/
cd $PROJECT_HOME/analysis/pca

for REF in Dog Wolf; do
    if [ $REF == "Dog" ]; then
	datadir=$DATA_HOME/dogRef
	curfai=$DOGFAI
    elif [ $REF == "Wolf" ]; then
	datadir=$DATA_HOME/wolfRef
	curfai=$WOLFFAI
    fi
    for MISSPCA in 1.0 0.95 0.9; do
	if [ ! -e $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.vcf ]; then
    	    bcftools view $datadir/allAlignedTo$REF.miss$MISSING.maf$MAF.minq30.mingq20.forPCA.bcf | vcftools --vcf - --max-missing $MISSPCA --recode --stdout > $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.vcf
	fi
	nsites=$(bcftools plugin counts $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.vcf | grep SNPs | cut -f2 -d:)
	nsamps=$(bcftools plugin counts $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.vcf | grep samples | cut -f2 -d:)
	if [ ! -e $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.geno ]; then
	    angsd -vcf-gl $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.vcf -fai $curfai -nInd $nsamps -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -out $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA && gunzip $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.geno.gz
	fi
	if [ ! -e $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.dosage ]; then
	    python $VCF2DOSAGE -v $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.vcf -o $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.dosage -m NA
	fi
	for MINMAF in 0.05 0.1 0.2; do 
	    ## NGSCovar used to generate pca plots
	    if [ ! -e $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MINMAF.minq30.mingq20.forPCA.covar ]; then
		$NGSCOVAR -probfile $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.geno -outfile allAlignedTo$REF.miss$MISSPCA.maf$MINMAF.minq30.mingq20.forPCA.normed.covar -nind $nsamps -nsites $nsites -call 0 -minmaf $MINMAF -norm 1
	    fi
	    ## Run eigen decomposition on the covar matrix using R
	    Rscript $PROJECT_HOME/code/pca.ngscovar.R allAlignedTo$REF.miss$MISSPCA.maf$MINMAF.minq30.mingq20.forPCA.normed.covar knownDatasets.samples.txt allAlignedTo$REF.miss$MISSPCA.maf$MINMAF.minq30.mingq20.ngscovar.pca.pdf
	    echo "========= MAF $MINMAF."
	done
	## Estimate the folded SFS from these files, using dosage files.
	echo "Done with MISS $MISSPCA."
    done
    ## Perform a pca using dosage and mds on the dosage files.
    Rscript $PROJECT_HOME/code/pca.dosage.R $datadir/allAlignedTo$REF.miss$MISSPCA.maf$MAF.minq30.mingq20.forPCA.dosage knownDatasets.samples.txt $MINMAF allAlignedTo$REF.miss$MISSPCA.maf$MINMAF.minq30.mingq20.dosage.pca.pdf
    echo "Done with $REF."
done

############################################
# Section to generate only matched reads & #
# do pca with only the matched reads. These#
# are the set of reads that are mapped uni-#
# quely in both species.                   #
############################################

cd $DATA_HOME
if [ ! -e matchedBams ]; then
    mkdir matchedBams
fi

#cd matchedBams
#$PROJECT_HOME/code/generateMatchedBams.sh

################################################
# Estimate SFS from the dosage file and angsd. #
################################################
