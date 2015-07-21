#! /usr/bin/env python

# File: plot_psmc_wolfdog.py
# Author: Shyam Gopalakrishnan
# Date: Jul 19 2015
# This script is used to plot the PSMC results 
# from alignment to dog and wolf. This requires 
# all the psmc results to be in one file and it
# also assumes that the order of samples in the
# file is alternating and matches. Looks like 
# sampl1dog, sample1wolf, sample2dog, sample2wolf
# and so on.

import argparse
import re
import sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy as np

def main():
    """This function is the main function that generates the plot.
    """
    psmcfile = open(args.psmcfile)
    labfile = open(args.labelfile)
    figure(1, figsize=(10,6))
    xmax = 0
    xmin = args.minx
    ymax = 0
    colors = ['b','g','r', 'c', 'm', 'y', 'k']
    linestyles=['-',':']
    curindex = 0
    while True:
        (times, popsizes) = processSingleSample(psmcfile, args.mut, args.gentime, args.round, args.window)
        if len(times) == 0:
            break
        if np.max(times) > xmax: xmax = np.max(times)
        if np.max(popsizes) > ymax: ymax = np.max(popsizes)
        label = labfile.readline()
        label = label.strip()
        print 'Using color', colors[int(curindex/2.0)%7], 'and ls', linestyles[curindex%2]
        step(times, popsizes, label=label, c=colors[int(curindex/2.0)%7], ls=linestyles[curindex%2])
        curindex += 1
    legend(loc=2, prop={'size':5})
    xlim([args.minx, xmax])
    ylim([1, args.maxy])
    ylabel('Population size')
    xlabel(r'Time ($\mu$='+str(args.mut)+', generation time='+str(args.gentime)+'yrs)')
    title('Effect of reference genome on PSMC results: Dog vs. Wolf')
    ax = gca()
    ax.set_xscale('log')
#    ax.set_yscale('log')
    savefig(args.out, bbox='tight')

def processSingleSample(infile, mut, gentime, chosenRound, window):
    """This function processes one part of the file
    to get the information for a single sample.
    The plotting itself is done in a separate function.
    """
    popsizes = []
    times = []
    rd = -1
    for line in infile:
        if re.match('RD', line) != None:
            rd = int(line.strip().split()[1])
        if rd != chosenRound: continue
        if re.match('TR', line) != None:
            theta = float(line.strip().split()[1])
            N0 = theta/(4*mut*window)
        if re.match('RS', line) != None:
            toks = [float(x) for x in line.strip().split()[2:4]]
            popsizes.append(toks[1]*N0)
            times.append(2*N0*gentime*toks[0])
        if re.match('PA', line) != None: break
    return (times, popsizes)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot psmc plots for wolf and dog reference project", version="0.1")
    parser.add_argument('-p', '--psmc', dest='psmcfile', type=str, required=True, help="Concatenated PSMC file")
    parser.add_argument('-l', '--label', dest='labelfile', type=str, required=True, help="Label file for legend")
    parser.add_argument('-m', '--mutation', dest='mut', type=float, help="Per generation per base mutation rate", default=1e-8)
    parser.add_argument('-g', '--gentime', dest='gentime', type=float, help="Generation time in years", default=3)
    parser.add_argument('-r', '--round', dest='round', type=int, help="Round to process", default=20)
    parser.add_argument('-w', '--window', dest='window', type=int, help="Window size while preparing data", default=100)
    parser.add_argument('--minX', dest='minx', type=float, help='Minimum x axis', default=1e4)
    parser.add_argument('--maxY', dest='maxy', type=float, help='Maximum population size', default=6e4)
    parser.add_argument('-o', '--output', dest='out', type=str, help="Output pdf file name", default="testout")
    args = parser.parse_args()
    main()
