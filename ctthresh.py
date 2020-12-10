#! /usr/bin/python
import sys
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import warnings
from rpy2.rinterface import RRuntimeWarning
from collections import Counter
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import os

warnings.filterwarnings("ignore", category=RRuntimeWarning)

def parse_user_input():
	"""
	Get and parse user input.

	"""

	parser = argparse.ArgumentParser()
	parser.add_argument('-f','--filter-file',required=True,help='Name of file containing filtered molecules (e.g. *exon_address.cfilt.txt or cellxgene file from kb)')
	parser.add_argument('-c','--cell-file',required=False,help='List of cell barcodes from kb.')
	parser.add_argument('-g','--gene-file',required=False,help='List of gene names from kb.')
	parser.add_argument('-r','--reference-file',required=True,help='Two-column file with gids and gene symbols')
	parser.add_argument('-p','--prefix',required=True,help='Prefix for output files including path')
	parser.add_argument('-l','--lower-thresh',required=True,help='Lower count threshold parameter for emptyDrops')
	return parser

# parse user input
parser = parse_user_input()
user_input = parser.parse_args()


cfilt_INFILE = user_input.filter_file 
cell_INFILE = user_input.cell_file
gene_INFILE = user_input.gene_file
ref_INFILE = user_input.reference_file
output_PREFIX = user_input.prefix
lower_THRESH = int(user_input.lower_thresh)
cbcind_OUTFILE = output_PREFIX+'.edrops.cbcs.txt'

bclist = [line.split()[0] for line in open(cell_INFILE)]
glist = [line.split()[0] for line in open(gene_INFILE)]
cbcs = {}
with open(cfilt_INFILE) as f:
	next(f)
	next(f)
	next(f)
	next(f) 
	for line in f:
		llist = line.split()
		cbc = bclist[int(llist[0])-1]
		for z in range(int(llist[2])):
			if cbc in cbcs.keys():
				cbcs[cbc].append(glist[int(llist[1])-1])
			else:
				cbcs[cbc] = [glist[int(llist[1])-1]]

for cbc in bclist:
	if cbc in cbcs.keys():
		if len(cbcs[cbc]) < lower_THRESH:
			del cbcs[cbc]

# get GID:gene symbol
gdict = {}
with open(ref_INFILE) as f:
	for line in f:
		llist = line.split()
		gdict[llist[0]] = llist[1]

# generate count matrix and summary statistics for kept cell barcodes
matrix = []
gcts = []
mcts = []
for cbc in cbcs.keys():
	cts = Counter(cbcs[cbc])
	matrix.append([cts[gid] if gid in cts.keys() else 0 for gid in gdict.keys()])
	gcts.append(len(set(cbcs[cbc])))
	mcts.append(len(cbcs[cbc]))
matrix = np.transpose(np.array(matrix))

# write count matrix to file
matrix_OUTFILE = output_PREFIX+'.matrix.txt'
with open(matrix_OUTFILE,'w') as g:
	for gid,vec in zip(gdict.keys(),matrix):
		gene = gdict[gid]
		st = gid+'\t'+gene+'\t'+'\t'.join([str(pt) for pt in vec])+'\n'
		g.write(st)

# write summary statistics (molecule and gene counts) to file
hist_OUTFILE = output_PREFIX+'.matrix.hist.txt'
with open(hist_OUTFILE,'w') as g:
	for cbc,mct,gct in zip(cbcs.keys(),mcts,gcts):
		st = cbc+'\t'+str(mct)+'\t'+str(gct)+'\n'
		g.write(st)


