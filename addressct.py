#! /usr/bin/python3.4
import sys

# The purpose of this program is to collapse read addresses based on UMI sequences without accounting for errors.

ct_gene_OUTFILE = sys.argv[1] # output for whole-gene alignment
ct_exon_OUTFILE = sys.argv[2] # output for exon-only alignment

address_dict = {}
for line in sys.stdin:	# file with read addresses - readid, cell barcode, umi, gene name, exon-only indicator
	llist = line.split()
	if len(llist) == 5:
		address = "\t".join(llist[1:4])
		case = int(llist[4])
		if case == 0: # ambiguous gene, unambiguous exon
			if address in address_dict.keys():
				address_dict[address][1]+=1
			else:
				address_dict[address] = [0,1]
		elif case == 1: # unambiguous gene, unambiguous exon
			if address in address_dict.keys():
				address_dict[address][0]+=1
				address_dict[address][1]+=1
			else:
				address_dict[address] = [1,1]
		elif case == 2: # unambiguous gene, no exon
			if address in address_dict.keys():
				address_dict[address][0]+=1
			else:
				address_dict[address] = [1,0]
					
output1 = open(ct_gene_OUTFILE,'w')
output2 = open(ct_exon_OUTFILE,'w')
for address in address_dict.keys():
	pt1 = address_dict[address][0]
	pt2 = address_dict[address][1]
	if pt1 > 0:
		output1.write('%(address)s\t%(pt1)d\n' % vars())
	if pt2 > 0:
		output2.write('%(address)s\t%(pt2)d\n' % vars())
output1.close()
output2.close()
