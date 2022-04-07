#! /usr/bin/python3.4
import sys
import gc
from collections import defaultdict
from bx.intervals.intersection import IntervalTree, Interval

# The purpose of this program is to extract reads that uniquely and strand-specifically align to genes in a GTF. The GTF must include a "gene" region category 
# for this program to work, and the alignment is assumed to have been conducted using STAR, although other aligners may be compatible.  The sam file should be
# piped into the program as stdin and the resulting address file should be piped out as stdout. It is highly recommended that the sam file be piped in using 
# "samtools view" (i.e. that the sam file actually be stored as a bam file) and that stdout be piped to gzip due to the likely size of these two files. It is also
# assumed that the sam file be sorted by genomic position (coordinate). This allows the annotation to be loaded into memory one chromosome at-a-time. Finally,
# it is assumed that the entire barcode sequence (e.g. cell/sample barcode + UMI) be included as a single string appended to the read ID in the sam file
# (separated from the read ID with a colon, e.g. machine_info:readid:barcode_sequence). The simplest way to achieve this is to append the barcode sequence to the
# read ID in the template fastq file prior to alignment. Note that genes are defined here to include introns and exons of all annotated isoforms.

# Updated 9/22/2016 to use interval trees for nucleotide lookup (get_interval_tree_for_chromosome function from Hanna)
# Updated 9/22/2016 to recognize multi-mappers with flag = 0, 16, 256, 272 (now counting primary mutli-mapped alignments)
# Updated 9/22/2016 to use new dictionary format for storing gtf (metadictionary containing a default dictionary for each strand from Hanna)
# Updated 9/23/2016 to incorporate simultaneous testing for whole-gene alignment and exon-only alignment
# Updated 6/14/2017 to accommodate 10x v2 libraries

def get_interval_tree_for_chromosome(gtf, chrm):
#   Returns interval trees of genes on the given chromosome for each strand,
#   arranged in a {strand: tree} dictionary
	tree_dict = {strand: IntervalTree() for strand in gtf.keys()}
	for strand, chrm_dict in gtf.items():
		for iv in chrm_dict[chrm]:
            # assumes iv of form (left, right, gene_id) as specified in get_gtf
			tree_dict[strand].insert_interval(Interval(*iv))	
	return tree_dict

gtf_INFILE = sys.argv[1] # preferable a Gencode-style gtf
technology = sys.argv[2] 

# For both the + and - strands, there is a default dictionary in the gtf metadictionary with chromosome names as keys and values
# that are lists of tuples (gene_id, left coordinate, right coordinate) indicate gene regions on that chromosome and strand.
genegtf = {'+':defaultdict(list),'-':defaultdict(list)} # whole-gene metadictionary keyed by '+' and '-' (indicating strand)
exongtf = {'+':defaultdict(list),'-':defaultdict(list)} # exon-only metadictionary keyed by '+' and '-' (indicating strand)
with open(gtf_INFILE) as g:
	for line in g:
		if line[0] != '#':
			llist = line.split()
			ch = llist[0]   # chromosome
			x1 = int(llist[3]) # left position
			x2 = int(llist[4]) # right position
			d = llist[6]       # strand (+ or -)
			gid = llist[llist.index('gene_id')+1][1:-2] # gene_id
			region = llist[2] # region type
			if region == 'gene':
				genegtf[d][ch].append((x1,x2+1,gid)) # load gene region into whole-gene metadictionary
			elif region == 'exon':
				exongtf[d][ch].append((x1,x2+1,gid)) # load exon region into exon-only metadictionary
				
current_ch = ''
store_mm_gene = {}
store_mm_exon = {}
genetree = ''
exontree = ''
if technology == 'DropSeqv1':
	cellbclen = 12
	umilen = 8
elif technology == '10xv2':
	cellbclen = 16
	umilen = 10
elif technology == 'DropSeqv2' or technology=='PearSeq':
	cellbclen = 12
	umilen = 9
elif technology == '10xv3':
	cellbclen = 16 
	umilen = 12

for line in sys.stdin:  # loop through coordinate-sorted sam file piped through stdin
	if line[0] != '@':
		llist = line.split()
		ch = llist[2] # chromosome
		flag = int(llist[1])  # alignment flag
		p1 = int(llist[3])  # left-most coordinate of alignment
		p2 = p1+len(llist[9])-1 # right-most coordinate of alignment
		idstring = llist[0].split(':')
		readid = ':'.join(idstring[3:7]) # parse read id
		cellbc = idstring[7][0:cellbclen]       # get cell barcode (first 12 or 16 nts)
		umi = idstring[7][cellbclen::] # get UMI (last 8 or 10 nts)
		hits = llist[11] # number of hits
		if ch != current_ch:  # if a new chromosome is encountered, load its nucleotides into memory
			del genetree
			gc.collect()
			del exontree
			gc.collect()
			genetree = get_interval_tree_for_chromosome(genegtf, ch)
			exontree = get_interval_tree_for_chromosome(exongtf, ch)
			current_ch = ch
		if flag == 0 or flag == 256:
			strand = '+'
		elif flag == 16 or flag == 272:
			strand = '-'
		gene_interval = genetree[strand].find(p1,p2) # try to get gene from + or - whole-gene tree
		exon_interval = exontree[strand].find(p1,p2) # try to get gene from + or - exon-only tree
		exon_gids = [iv.value for iv in exon_interval] # how many genes are associated with the exon intervals?
		n_exon_gids = len(set(exon_gids)) # need to count because multiple exons per gene
		n_gene_intervals = len(gene_interval) # only one "gene per gene", don't need to count gene ids for whole-gene
		if n_gene_intervals > 1:
			if n_exon_gids == 1:	# case 0: whole-gene is ambiguous but exon-only is unambiguous
				gid = exon_gids[0]
				if hits == 'NH:i:1':
					if len(umi) == umilen:
						sys.stdout.write('%(readid)s\t%(cellbc)s\t%(umi)s\t%(gid)s\t0\n' % vars()) # stdout piped to address file (preferably through gzip), case 0 
				else: # handle multi-mappers- will keep multimappers if there's only one that maps strand-specifically and within a gene
					if readid not in store_mm_exon: 
						store_mm_exon[readid] = gid,cellbc,umi,0
					else:
						store_mm_exon[readid] = -1 # if readid is already in store_mm, then there must be more than one strand-specific alignment, so discard
		elif n_gene_intervals == 1:     
			if n_exon_gids == 1:	# case 1: whole-gene is unambiguous, therefore exon-only is unambiguous if it exists
				gid = exon_gids[0]
				if hits == 'NH:i:1':
					if len(umi) == umilen:
						sys.stdout.write('%(readid)s\t%(cellbc)s\t%(umi)s\t%(gid)s\t1\n' % vars()) # stdout piped to address file (preferably through gzip), case 1
				else:
					if readid not in store_mm_exon:
						store_mm_exon[readid] = gid,cellbc,umi,1
					else:
						store_mm_exon[readid] = -1
					if readid not in store_mm_gene:
						store_mm_gene[readid] = gid,cellbc,umi,1
					else:
						store_mm_gene[readid] = -1
			elif n_exon_gids == 0:	# case 2: whole gene is unambiguous and exon-only does not exist
				gid = gene_interval[0].value
				if hits == 'NH:i:1':
					if len(umi) == umilen:
						sys.stdout.write('%(readid)s\t%(cellbc)s\t%(umi)s\t%(gid)s\t2\n' % vars()) # stdout piped to address file (preferably through gzip), case 2
				else:
					if readid not in store_mm_gene:
						store_mm_gene[readid] = gid,cellbc,umi,2
					else:
						store_mm_gene[readid] = -1
				
			
for readid in store_mm_gene.keys(): # write all multi-mappers that are actually unique
	if store_mm_gene[readid] != -1:
		gid = store_mm_gene[readid][0]
		cellbc = store_mm_gene[readid][1]
		umi = store_mm_gene[readid][2]
		case = store_mm_gene[readid][3]
		if len(umi) == umilen:
			sys.stdout.write('%(readid)s\t%(cellbc)s\t%(umi)s\t%(gid)s\t%(case)d\n' % vars())
for readid in store_mm_exon.keys(): # write all multi-mappers that are actually unique
	if store_mm_exon[readid] != -1:
		gid = store_mm_exon[readid][0]
		cellbc = store_mm_exon[readid][1]
		umi = store_mm_exon[readid][2]
		case = store_mm_exon[readid][3]
		if len(umi) == umilen:
			sys.stdout.write('%(readid)s\t%(cellbc)s\t%(umi)s\t%(gid)s\t%(case)d\n' % vars())

