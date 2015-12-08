#! /bin/bash

# File: generateMatchedBams.sh
# Author: SG
# Description: Generates bams that contain only
# reads that aligned uniquely are present in both bams.
# Uses bedtools intersect to accomplish this.

module load bedtools/
