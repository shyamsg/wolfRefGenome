#! /usr/bin/env python

#####################################
# File: computeAvgDepthFromStats.py #
# Author: Shyam Gopalakrishnan      #
# Date: 06/10/2015                  #
# Description: Computes the mean    #
# depth of the sample given the stat#
# file for that sample generated    #
# using bcftools.                   #
#####################################

import numpy as np
import sys
import argparse as ap

if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Depth computer from bcftools stats file.")
    parser.add_argument("-s", "--statfile", type=str, metavar="StatFile", dest="statfile", help="Statistics file for sample")
    args = parser.parse_args()
    statfile = args.statfile
else:
    statfile = "test.stats"

cumulDepth = 0.0
numSites = 0
for line in infile:
    line = line.strip().split()
    if line[0] == "DP":
        curSites = int(line[5])
        numSites += curSites
        if line[2][0:1] == ">":
            cumulDepth += curSites*(float(line[2][1:])+1)
        else:
            cumulDepth += curSites*float(line[2])
            
print statfile[0:-6], np.round(cumulDepth/numSites, 3)
