#! /usr/bin/python3.4
import sys
import operator
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from numpy import mean

# The purpose of this program is to generate some figures containing
# summary statistics about a DropSeq experiment.


pfilt_INFILE = sys.argv[1] # the output of DropSeqPipeline5_pfilter.py
pdf_OUTFILE = sys.argv[2]
chist_OUTFILE = sys.argv[3]

print('Loading filtered data...')
bc_mdict = {}
bc_gdict = {}
with open(pfilt_INFILE) as f:
	for line in f:
		llist = line.split()
		bc = (llist[0])
		if bc in bc_mdict:
			bc_mdict[bc]+=1 # Number of molecules for each cell barcode
			bc_gdict[bc].append(llist[2]) # list of genes for each cell barcode
		else:	
			bc_mdict[bc]=1
			bc_gdict[bc] = [llist[2]]

for bc in bc_gdict.keys():
	bc_gdict[bc] = len(list(set(bc_gdict[bc]))) # collapse to a list of unique genes

Nbc = len(bc_gdict.keys())
print('Found %(Nbc)d cell barcodes...' % vars())
print('Sorting cell barcodes...')
sorted_mdict = sorted(bc_mdict.items(), key=operator.itemgetter(1), reverse=True) # sort cell barcodes by number of molecules

msum = 0
for i in range(0,len(sorted_mdict)):
	msum+=sorted_mdict[i][1]
msum=float(msum) # total molecules in the data set

print('Generating cumulative histogram and summary statistics...')
output = open(chist_OUTFILE,'w')
mchist = []
pt3 = 0.0
molecules = []
avgmolecules = []
genes = []
avggenes = []
m=0
g=0
j=1
for pt in sorted_mdict: # compute cumulative histogram for molecules per cell barcode
	pt1 = pt[0]
	pt2 = pt[1]
	pt3 += float(pt2)/msum
	mchist.append(pt3)
	molecules.append(pt2)
	pt4=bc_gdict[pt1]
	genes.append(pt4)
	m+=pt2
	g+=pt4
	avgmolecules.append(float(m)/float(j))
	avggenes.append(float(g)/float(j))
	output.write('%(pt1)s\t%(pt2)d\t%(pt3)f\n' % vars())
	j=j+1
output.close()

print('Generating pdf...')
with PdfPages(pdf_OUTFILE) as pdf:
	plt.plot(mchist)
	plt.xlabel('Cell Barcode')
	plt.ylabel('Cumulative Histogram for Molecules')
	pdf.savefig()
	plt.close()
	plt.plot(mchist)
	plt.xlabel('Cell Barcode')
	plt.ylabel('Cumulative Histogram for Molecules')
	plt.xscale('log')
	pdf.savefig()
	plt.close()	
	plt.plot(molecules)
	plt.xlabel('Cell Barcode')
	plt.ylabel('Molecules')
	plt.xscale('log')
	plt.yscale('log')
	pdf.savefig()
	plt.close()
	plt.plot(avgmolecules)
	plt.xlabel('Number of Cells Included')
	plt.ylabel('Average Number of Molecules per Cell')
	plt.xscale('log')
	plt.yscale('log')
	pdf.savefig()
	plt.close()
	plt.plot(avggenes)
	plt.xlabel('Number of Cells Included')
	plt.ylabel('Average Number of Genes per Cell')
	plt.xscale('log')
	plt.yscale('log')
	pdf.savefig()
	plt.close()
