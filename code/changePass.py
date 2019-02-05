#! /bin/bash

import gzip as gz
import sys

f=gz.open(sys.argv[1])
o=open(sys.argv[2], "w")

for line in f:
    if line[0] != "#":
        toks = line.strip().split()
        for i in [2,3,4,5,7,8,9]:
            toks[i] = toks[i].replace('PASS', '.')
        line = "\t".join(toks)+'\n'
    o.write(line)
f.close()
o.close()
