#! /usr/bin/env python

# File: plot_sampleStats.py
# Author: Shyam Gopalakrishnan
# Date: 22nd July 2015
# This script is used to process and plot depth
# and geno quality statistics from a single or a 
# multi sample vcf. 

import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy as np
import sys
import argparse
import gzip

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot sample statisitcs from vcf", version="0.1")
    parser.add_argument('-v', '--vcf', metavar='VCFFile', dest='vcf', help='Input vcf file (or gzvcf file)', )
