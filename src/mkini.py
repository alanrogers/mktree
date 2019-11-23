#!/usr/bin/python

# mkini reads an input file with one column per epoch and writes one
# in the format of the file pophist.ini. The input of mkini should have
# 4 columns:
#
# col 1: haploid size of each subdivision
# col 2: mean number of immigrants per generation into each subdivision
# col 3: time in generations of the endpoint (closest to root) of the
#        current epoch.   
# col 4: number of subdivisions

import sys

# mutation rate per sequence per generation. In other words, the sum across
# sites of the rate per site per generation.
mut = 1e-3

if len(sys.argv) != 2:
    print "usage: mkini.py input_file"
    sys.exit(1)
try:
    f = open(sys.argv[1],"r")
except:
    print "Couldn't open file \"%s\"" % sys.argv[1]
    sys.exit(1)

thusfar = 0.0

print "%%%7s %8s %8s %8s" % ("theta", "mn", "tau", "npops")
for line in f:
    line=line.strip().split()
    theta = float(line[0]) * 2.0 * mut
    mn = float(line[1])
    t = float(line[2])  # should handle "inf" gracefully
    tau = (t - thusfar) * 2.0 * mut
    thusfar = t
    npops = int(line[3])
    print "%8.4f %8.4f %8.4f %d" % (theta, mn, tau, npops)
    
