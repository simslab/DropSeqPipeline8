#! /usr/bin/python3.4
import sys

# The purpose of this program is to collapse read addresses based on UMI sequences without accounting for errors.
# The program is used for a second round of UMI collapsing that takes place after cell barcode correction. The cell
# barcode correction alters the cell barcode and UMI sequence, and so a second round of collapse is required.

ct_OUTFILE = sys.argv[1] # output

address_dict = {}
for line in sys.stdin:	# file with molecule addresses - cell barcode, umi, gene name, read count
	llist = line.split()
	if len(llist) == 4:
		address = "\t".join(llist[0:3])
		if address in address_dict.keys():
			address_dict[address] += int(llist[3])    # count the number of times a molecule address occurs
		else:
			address_dict[address] = int(llist[3])

output = open(ct_OUTFILE,'w')
for address in address_dict.keys():
	pt = address_dict[address]
	output.write('%(address)s\t%(pt)d\n' % vars())
output.close()

