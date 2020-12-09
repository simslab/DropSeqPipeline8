#! /usr/bin/python3.4
import sys
import gzip

# The first purpose of this program is to extract the cell barcode and UMI from the read 1 fastq file of a single cell RNA-Seq experiment
# and append it to the read ID of the corresponding line in the read 2 fastq file.  The second purpose is to clip the poly(A)-tail from 
# 3'-end of read 2 and discard any sequences that are too short for alignment after clipping. The function below is designed to be compiled
# with cython. The gzipped read 1 fastq is the first argument and the second input should be the read 2 fastq piped in from stdin (e.g. with
# zcat).  The resulting read 2 fastq file is written to stdout and should be compressed externally.

def clipper(fastq1_INFILE,technology,fastq2_INFILE): # gzip compressed input and output
	bc_seqs = [] # list of sequences containing the 12-nt cell barcode followed by the 8-nt UMI
	bclen = 0
	if technology == '10xv2':
		bclen = 26 # 10x v2 barcode region is first 26 nt of read 1
	elif technology == 'DropSeqv1':
		bclen = 20 # drop-seq barcode region is first 20 nt of read 1
	elif technology == 'DropSeqv2':
		bclen = 21 # new drop-seq barcode region is first 21 nt of read 1
	elif technology == '10xv3':
		bclen = 28 # 10x v3 barcode region is first 28 nt of read 1
	i = 0
	with gzip.open(fastq1_INFILE,'rb') as f:
		for line in f:
			if i == 0:
				i = 1
			elif i == 1:
				bc_seqs.append(line[0:bclen].decode()) # first bclen nts of read 1 contain cell barcode and UMI
				i = 2
			elif i == 2:
				i = 3
			elif i == 3:
				i = 0
	go = 0
	a = 'AAAAAAAA' # poly(A) tail
	j = 0
	i = 0
	for line in fastq2_INFILE: # from stdin
		if i == 0:
			if bc_seqs[j][-1] == '\n':
				bc_seqs[j] = bc_seqs[j][:-1]+'N'
			line1 = line.split()[0]+':'+bc_seqs[j]+'\n' # add 20 nt barcode sequence to read ID separated by colon
			j+=1
			go = 0
			i = 1
		elif i == 1:
			if a in line: # if poly(A) tail is in read 2
				x  = line.find(a) # find its position
				if x > 23:  # if clipped read is sufficiently long
					line2 = line[0:x]+'\n' # clip it and keep
					go = 1
			else: # also keep it if there's no poly(A) tail
				x = -1
				line2 = line 
				go = 1
			i = 2
		elif i == 2:
			line3 = line # comment line
			i = 3
		elif i == 3: # q-score line
			if go == 1:
				if x != -1: 
					newlines = line1+line2+line3+line[0:x]+'\n'
				else:
					newlines = line1+line2+line3+line
				sys.stdout.write(newlines)
			i = 0
	return 0

# Enumerate HD=1 sequences for a given barcode sequence
def enumerate_bcs(bc):
	bcs = []
	N = len(bc)
	for i in range(N):
		bc1=''
		bc2=''
		bc3=''
		bc4=''
		bc5=''
		for j in range(N):
			if i==j:
				bc1+='A'
				bc2+='G'
				bc3+='C'
				bc4+='T'
				bc5+='N'
			else:
				bc1+=bc[j]
				bc2+=bc[j]
				bc3+=bc[j]
				bc4+=bc[j]
				bc5+=bc[j]
		bcs.append(bc1)
		bcs.append(bc2)
		bcs.append(bc3)
		bcs.append(bc4) 
		bcs.append(bc5)
	bcs = list(set(bcs))
	return bcs


def pearclipper(fastq1_INFILE,fastq2_INFILE,pearaddress_OUTFILE,citeseq_INFILE): # gzip compressed input and output
	bc_seqs = [] 
	ext_citeseq = {}
	with open(citeseq_INFILE) as f:
		for line in f:
			llist = line.split()
			bc = llist[1]
			bcs = enumerate_bcs(bc)
			for bc2 in bcs:
				ext_citeseq[bc2] = llist[0]
	bclen = 21
	cbclen = 12
	cslen = len(list(ext_citeseq.keys())[0])
	i = 0
	with gzip.open(fastq1_INFILE,'rb') as f:
		for line in f:
			if i == 0:
				i = 1
			elif i == 1:
				bc_seqs.append(line[0:bclen].decode()) # first bclen nts of read 1 contain cell barcode and UMI
				i = 2
			elif i == 2:
				i = 3
			elif i == 3:
				i = 0
	j = 0
	i = 0
	with open(pearaddress_OUTFILE,'w') as g:
		for line in fastq2_INFILE: # from stdin
			if i == 0:
				if bc_seqs[j][-1] == '\n':
					bc_seqs[j] = bc_seqs[j][:-1]+'N'
				llist = line.split()
				readid = ':'.join(llist[0].split(':')[3:7])
				go = 0
				i = 1
			elif i == 1:
				csbc = line[0:cslen]
				if csbc in ext_citeseq.keys():
					if len(bc_seqs[j]) == 21:
						cbc = bc_seqs[j][0:12]
						umi = bc_seqs[j][12::]
						feature = ext_citeseq[csbc]
						g.write('%(readid)s\t%(cbc)s\t%(umi)s\t%(feature)s\t0\n' % vars())
				i = 2
				j+=1
			elif i == 2:
				i = 3
			elif i == 3: # q-score line
				i = 0
	return 0

