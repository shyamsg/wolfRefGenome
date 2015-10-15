#! /usr/bin/env python

# File: filterSampleLevel.py
# Date: 22nd July 2015
# Author: Shyam Gopalakrishnan
# This script removes the null indels, indels with
# . as the alternate allele, meaning no alt allele.
# This script filters sites based on a couple of 
# criteria -
# 1. distance from indel
# 2. distance from other snvs
# 3. genotype qualities
# The repeat mask removal and triallelic 
# removal is done after this step.
# Once the samples are all merged into a single
# vcf, you can determine if the site is triallelic
# or not. Also, using vcftools to get only a subset
# of regions is easy, so repeat mask removal is also
# done after this script. 

import re
import sys
import argparse
import gzip

def filterSample(vcf, out, qual, nQual, vardist, mind, maxd):
    """The function filters the sample level file to remove 
    null indels, and variants that fail certain quality
    constraints, such as distance to closest SNP or indel 
    minimum and maximum depth of coverage. 
    """
    if vcf[-2:] == "gz":
        f = gzip.open(vcf)
    else:
        f = open(vcf)

    if out == "stdout":
        outfile = sys.stdout
    else:
        outfile = open(out, "w")

    while True:
        line = f.readline()
        if line[0] == "#":
            outfile.write(line)
        if line.strip().split()[0] == '#CHROM':
            break
        
    cnt = 0
    lines = []
    linepos = []
    linetype = []
    linequal = []
    linegeno = []
    curChr = ''
#    print "0 lines done"
    for line in f:
        cnt += 1 
        if cnt % 1000000 == 0:
            print "\r", cnt, "lines done"
        toks = line.strip().split()
        if toks[0] != curChr: # next chromosome
            curChr = toks[0]
            for l2 in lines:
                outfile.write(l2)
            del lines[:]
            del linepos[:]
            del linequal[:]
            del linegeno[:]
            del linetype[:]
        curtype = 0
        if re.match("INDEL", toks[7]) != None:
            if toks[4] == ".":
                # Skip if null indel
#                print line, "null indel"
                continue
            else: 
                curtype = 2 #INDEL
        elif toks[4] != '.':
            curtype = 1 # SNP
        curpos = int(toks[1])
        curdepth = int([x.split('=')[1] for x in toks[7].split(';') if x[0:2] == 'DP'][0])
        # site fails depth criteria
        if curdepth < mind or curdepth > maxd:
#            print curpos, curdepth, "depth filter"
            continue
        # if curline is not variant, add it and continue
        if curtype == 0:
            lines.append(line)
            linequal.append(0)
            linetype.append(curtype)
            linepos.append(curpos)
            linegeno.append('')
            continue
        toks2 = toks[9].split(':')
        curqual = int(toks2[2])
        curgeno = toks2[0]
        # if curqual is low, set getnotype to missing
        if curqual < qual:
            line = line.replace(curgeno, './.', 1)
        # Since curline is variant, go over previous stuff and output the stuff 
        # outside range and delete them later
        delindexes = []
        if curtype == 1: # snp
            for index in xrange(len(lines)):
                # sites that are not influenced by new site
                if (linepos[index]+vardist < curpos):
                    outfile.write(lines[index])
#                    print "snp write", linepos[index]
                    delindexes.append(index)
                    continue
                # Sites that are influenced by current site
                # and thus influence the current site
                # invariant has no effect, so only looking at snps and indels
                if linequal[index] < nQual:
                    continue
                if linetype[index] == 1: # SNP
                    lines[index] = lines[index].replace(linegeno[index], './.', 1)
                    line = line.replace(curgeno, './.', 1)
                elif linetype[index] == 2: # indel
                    line = line.replace(curgeno, './.', 1)
        else: # curtype has to be INDEL, since non-variants don't make it this far
            for index in xrange(len(lines)):
                # sites that are not influenced by new site
                if (linepos[index]+vardist < curpos):
                    outfile.write(lines[index])
#                    print "indel write", linepos[index]
                    delindexes.append(index)
                    continue
                # Sites that are influenced by current site
                # and thus influence the current site
                # invariant has no effect, so only looking at snps and indels
                if linequal[index] < nQual:
                    continue
                if linetype[index] == 1: # SNP
                    lines[index] = lines[index].replace(linegeno[index], './.', 1)
        for index in reversed(delindexes):
            del lines[index]
            del linepos[index]
            del linequal[index]
            del linegeno[index]
            del linetype[index]
        del delindexes[:]
        # add current line to lists
        lines.append(line)
        linepos.append(curpos)
        linetype.append(curtype)
        linegeno.append(curgeno)
        linequal.append(curqual)
    # end of file cleanup
    for index in xrange(len(lines)):
        outfile.write(lines[index])
        print "final write", linepos[index]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove null indels from vcf")
    parser.add_argument("-i", "--invcf", metavar="vcffile", dest="vcf", required=True, 
                        type=str, help="Input VCF file (gzipped or not)")
    parser.add_argument("-o", "--out", metavar="outfile", dest="out", required=True,
                        type=str, help="Output VCF file (unzipped) OR stdout")
    parser.add_argument("-q", "--qual", metavar="quality", dest="qual", required=False, 
                        type=int, default=20, help="Minimum genotype quality of variant (both SNP and indel)")
    parser.add_argument("-n", "--neighborQual", metavar="nQual", dest="nqual", required=False, 
                        type=int, default=10, help="Minimum quality of SNP/indel neighbor (to filter by distance)")
    parser.add_argument("-d", "--minDepth", metavar="minDepth", dest="mindepth", required=False, 
                        type=float, default=5, help="Minimum depth to retain variants")
    parser.add_argument("-D", "--maxDepth", metavar="maxDepth", dest="maxdepth", required=False, 
                        type=float, default=100, help="Maximum depth to retain variants")
    parser.add_argument("-v", "--varDist", metavar="varDistance", dest="vardist", 
                        required=False, type=int, default=5, help="Minimum distance from snp/indel")
    args = parser.parse_args()
    if args.qual < args.nqual:
        print "The quality threshold is less than neighboring quality threshold."
        print "This is not sensible, so setting neighboring quality to quality threshold."
        args.nqual = args.qual
    filterSample(args.vcf, args.out, args.qual, args.nqual, args.vardist, args.mindepth, args.maxdepth)
