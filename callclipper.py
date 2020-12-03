#! /usr/bin/python3.4
import sys
from DropSeqPipeline7_clipper import clipper,pearclipper

technology = sys.argv[2]
if technology == 'PearSeq':
    pearclipper(sys.argv[1],sys.stdin,sys.argv[3],sys.argv[4])
else:
    clipper(sys.argv[1],technology,sys.stdin) # call cython-compiled poly(A) clipper function- first argument is read 1 fastq (gzipped) and second argument is read 2 fastq piped from stdin
