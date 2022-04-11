#! /usr/bin/python
import os
import argparse
from random import random

def parse_user_input():
	"""
	Get and parse user input.
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument('-o','--outdir',required=True,help='Path to output directory.')
	parser.add_argument('-p','--prefix',required=True,help='Prefix for output files.')
	parser.add_argument('-r1','--read1fastq',required=True,help='Path to read 1 fastqs separated by commas.')
	parser.add_argument('-r2','--read2fastq',required=True,help='Path to read 2 fastqs searated by commas.')
	parser.add_argument('-r','--refdir',required=False,help='Path to reference genome directory.')
	parser.add_argument('-g','--gtf',required=False,help='Path to gtf.')
	parser.add_argument('-v','--overhang',type=int,required=False,help='sjdbOverhang parameter for STAR.')
	parser.add_argument('-t','--threads',type=int,required=True,help='Number of threads for STAR and samtools.')
	parser.add_argument('-s','--s3path',required=True,help='s3 bucket for data storage.')
	parser.add_argument('-x','--technology',required=True,choices=['10xv2','10xv3','DropSeqv1','DropSeqv2','PearSeq','CiteSeq3v3','CiteSeq5v2','CiteSeqTSB'],help='Technology (determines barcoding scheme).') 
	parser.add_argument('-b','--citeseq-barcodes',required=False,help='Path to 2-column file with cite-seq features and barcodes. Required if technology is PearSeq.')
	parser.add_argument('-e','--exon-only',action='store_true',help='Option to skip whole-gene body processing and only process exon-only reads.')
	parser.add_argument('-a','--address-only',action='store_true',help='Option to stop after generating an address file (e.g. for index-swap correction).')
	parser.add_argument('-ps','--post-swap',action='store_true',help='Option to start pipeline from a pre-generated address file (e.g. after index-swap correction).')
	return parser


parser = parse_user_input()
ui = parser.parse_args()

outdir = ui.outdir
prefix = ui.prefix
r1fastq_INFILES = ui.read1fastq.split(',')
r2fastq_INFILES = ui.read2fastq.split(',')
r2fastqclip_INFILE = outdir+'/'+prefix+'_R2.clip.fastq.gz'
technology = ui.technology
s3bucket = ui.s3path
if not ui.post_swap:
	if technology == 'PearSeq' or technology in ['CiteSeq3v3','CiteSeq5v2','CiteSeqTSB']:
		if not ui.citeseq_barcodes:
			print('Error: --citeseq-barcodes path required for PearSeq')
			exit()
		else:
			citeseq_INFILE = ui.citeseq_barcodes
			pearaddress_INFILE = outdir+'/'+prefix+'.pear_address.txt'
			i=0
			for r1fastq_INFILE,r2fastq_INFILE in zip(r1fastq_INFILES,r2fastq_INFILES):
				tmp = outdir+'/'+prefix+'.pear_address.'+str(i)+'.txt'
				cmd = 'zcat %(r2fastq_INFILE)s | python callclipper.py %(r1fastq_INFILE)s %(technology)s %(tmp)s %(citeseq_INFILE)s' % vars()
				print(cmd)
				os.system(cmd)
				i+=1
			cmd = 'cat %(outdir)s/%(prefix)s*.pear_address.*.txt | gzip > %(pearaddress_INFILE)s.gz' % vars()
			print(cmd)
			os.system(cmd)
			cmd = 'rm %(outdir)s/%(prefix)s*.pear_address.*.txt' % vars()
			print(cmd)
			os.system(cmd)
			if ui.address_only:
				print('Stopping after address file generation...')
				exit()
	else:
		i=0
		for r1fastq_INFILE,r2fastq_INFILE in zip(r1fastq_INFILES,r2fastq_INFILES):
			tmp = outdir+'/'+prefix+'_'+str(i)+'_R2.clip.fastq.gz'
			cmd = 'zcat %(r2fastq_INFILE)s | python callclipper.py %(r1fastq_INFILE)s %(technology)s | gzip > %(tmp)s' % vars()	 
			print(cmd)
			os.system(cmd)
			i+=1
		cmd = 'cat %(outdir)s/%(prefix)s_*_R2.clip.fastq.gz > %(r2fastqclip_INFILE)s' % vars()
		print(cmd)
		os.system(cmd)
#		cmd = 'rm %(outdir)s/%(prefix)s_*_R2.clip.fastq.gz' % vars()
#		print(cmd)
#		os.system(cmd)
		
		

# copy fastqs to S3
	for r1fastq_INFILE,r2fastq_INFILE in zip(r1fastq_INFILES,r2fastq_INFILES):
		cmd = 'aws s3 cp %(r1fastq_INFILE)s s3://%(s3bucket)s/%(r1fastq_INFILE)s' % vars()
		print(cmd)
		os.system(cmd)
		cmd = 'aws s3 cp %(r2fastq_INFILE)s s3://%(s3bucket)s/%(r2fastq_INFILE)s' % vars()
		print(cmd)
		os.system(cmd)

# remove fastqs
		cmd = 'rm %(r1fastq_INFILE)s %(r2fastq_INFILE)s' % vars()
		print(cmd)
		os.system(cmd)

# align clipped read 2 fastq to the genome/transcriptome annotation with 2-pass STAR
if technology not in ['PearSeq','CiteSeq5v2','CiteSeq3v3','CiteSeqTSB']:
	if not ui.post_swap:
		bam_OUTFILE = outdir+'/'+prefix+'.Aligned.out.bam'
		refdir = ui.refdir
		gtf_INFILE = ui.gtf
		threads = ui.threads
		overhang = ui.overhang
		cmd = 'python star.py %(outdir)s %(prefix)s %(r2fastqclip_INFILE)s %(refdir)s %(gtf_INFILE)s %(threads)d %(overhang)d' % vars()
		print(cmd)
		os.system(cmd)

		# sort bam file by coordinate
		sortedbam_OUTFILE = outdir+'/'+prefix+'.Aligned.out.sorted.bam'
		cmd = 'samtools sort -@ %(threads)d %(bam_OUTFILE)s -o %(sortedbam_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		# remove unsorted bam file
		cmd = 'rm %(bam_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		# extract read addresses that uniquely and strand-specifically align to whole genes (including introns)
		address_OUTFILE = outdir+'/'+prefix+'.address.txt.gz'
		cmd = 'samtools view %(sortedbam_OUTFILE)s | python address.py %(gtf_INFILE)s %(technology)s | gzip > %(address_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		# copy bam file to S3
		cmd = 'aws s3 sync %(outdir)s s3://%(s3bucket)s/%(outdir)s' % vars()
		print(cmd)
		os.system(cmd)

		# remove bam file
		cmd = 'rm %(sortedbam_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		if ui.address_only:
			print('Stopping after address file generation...')
			exit()

	else:
		address_OUTFILE = outdir+'/'+prefix+'.address.swap.decode.txt.gz'

	# count the number of times each read address occurs for whole genes (UMI collapse without error correction)
	geneaddresscts_OUTFILE = outdir+'/'+prefix+'.gene_address.cts.txt'
	exonaddresscts_OUTFILE = outdir+'/'+prefix+'.exon_address.cts.txt'
	cmd = 'zcat %(address_OUTFILE)s | python addressct.py %(geneaddresscts_OUTFILE)s %(exonaddresscts_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

	# collapse cell barcodes that arise from incomplete extension (truncated cell barcodes) and re-write their UMIs for whole genes
	if not ui.exon_only:
		geneaddressccts_OUTFILE = outdir+'/'+prefix+'.gene_address.ccts.txt'
		gene_trackbcs1_OUTFILE = outdir+'/'+prefix+'.gene_address.trackbcs1.txt'
		cmd = 'python collapse.py %(geneaddresscts_OUTFILE)s %(geneaddressccts_OUTFILE)s %(gene_trackbcs1_OUTFILE)s %(technology)s' % vars()
		print(cmd)
		os.system(cmd)

	# re-count the number of times each read address occurs for whole genes (UMI collapse without error correction)
		geneaddresscccts_OUTFILE = outdir+'/'+prefix+'.gene_address.cccts.txt'
		cmd = 'cat %(geneaddressccts_OUTFILE)s | python addressct2.py %(geneaddresscccts_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		# apply H=1 filter to eliminate UMIs that arise from sequencing error for whole genes
		genefiltcts_OUTFILE = outdir+'/'+prefix+'.gene_address.filt.txt'
		cmd = 'cat %(geneaddresscccts_OUTFILE)s | python callfilter.py %(genefiltcts_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		# clean up temporary files
		cmd = 'rm %(geneaddresscts_OUTFILE)s %(geneaddressccts_OUTFILE)s %(geneaddresscccts_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		# apply overlap-thresholded H=1 filter to collapse cell barcodes resulting from sequencing error (require at least 20 molecules and 75% molcular overlap for collapse)
		# for whole genes
		geneaddresscts2_OUTFILE = outdir+'/'+prefix+'.gene_address.cts2.txt'
		gene_trackbcs_OUTFILE = outdir+'/'+prefix+'.gene_address.trackbcs.txt'
		cmd = 'python collapse2.py %(genefiltcts_OUTFILE)s %(geneaddresscts2_OUTFILE)s 20 0.75 %(gene_trackbcs1_OUTFILE)s %(gene_trackbcs_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		cmd = 'rm %(gene_trackbcs1_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		# re-count the number of times each read address occurs for whole genes (UMI collapse without error correction)
		geneaddressccts2_OUTFILE = outdir+'/'+prefix+'.gene_address.ccts2.txt'
		cmd = 'cat %(geneaddresscts2_OUTFILE)s | python addressct2.py %(geneaddressccts2_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		# apply H=1 filter to eliminate UMIs that arise from sequencing error for whole genes
		genecfilt_OUTFILE = outdir+'/'+prefix+'.gene_address.cfilt.txt'
		cmd = 'cat %(geneaddressccts2_OUTFILE)s | python callfilter.py %(genecfilt_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		# clean up temporary files
		cmd = 'rm %(geneaddresscts2_OUTFILE)s %(geneaddressccts2_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

		# compute summary statistics on UMI-filtered data for whole genes (this time after cell barcode collapse)
		genecfiltpdf_OUTFILE = outdir+'/'+prefix+'.gene_address.cfilt.pdf'
		genecfiltchist_OUTFILE = outdir+'/'+prefix+'.gene_address.cfilt.chist.txt'
		cmd = 'python dsstats.py %(genecfilt_OUTFILE)s %(genecfiltpdf_OUTFILE)s %(genecfiltchist_OUTFILE)s' % vars()
		print(cmd)
		os.system(cmd)

	# collapse cell barcodes that arise from incomplete extension (truncated cell barcodes) and re-write their UMIs for exons
	exonaddressccts_OUTFILE = outdir+'/'+prefix+'.exon_address.ccts.txt'
	exon_trackbcs1_OUTFILE = outdir+'/'+prefix+'.exon_address.trackbcs1.txt'
	cmd = 'python collapse.py %(exonaddresscts_OUTFILE)s %(exonaddressccts_OUTFILE)s %(exon_trackbcs1_OUTFILE)s %(technology)s' % vars()
	print(cmd)
	os.system(cmd)

	# re-count the number of times each read address occurs for exons (UMI collapse without error correction)
	exonaddresscccts_OUTFILE = outdir+'/'+prefix+'.exon_address.cccts.txt'
	cmd = 'cat %(exonaddressccts_OUTFILE)s | python addressct2.py %(exonaddresscccts_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

	# apply H=1 filter to eliminate UMIs that arise from sequencing error for exons
	exonfiltcts_OUTFILE = outdir+'/'+prefix+'.exon_address.filt.txt'
	cmd = 'cat %(exonaddresscccts_OUTFILE)s | python callfilter.py %(exonfiltcts_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

	# clean up temporary files
	cmd = 'rm %(exonaddresscts_OUTFILE)s %(exonaddressccts_OUTFILE)s %(exonaddresscccts_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

	# apply overlap-thresholded H=1 filter to collapse cell barcodes resulting from sequencing error (require at least 20 molecules and 75% molcular overlap for collapse)
	# for exons
	exonaddresscts2_OUTFILE = outdir+'/'+prefix+'.exon_address.cts2.txt'
	exon_trackbcs_OUTFILE = outdir+'/'+prefix+'.exon_address.trackbcs.txt'
	cmd = 'python collapse2.py %(exonfiltcts_OUTFILE)s %(exonaddresscts2_OUTFILE)s 20 0.75 %(exon_trackbcs1_OUTFILE)s %(exon_trackbcs_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

	cmd = 'rm %(exon_trackbcs1_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

	# re-count the number of times each read address occurs for exons (UMI collapse without error correction)
	exonaddressccts2_OUTFILE = outdir+'/'+prefix+'.exon_address.ccts2.txt'
	cmd = 'cat %(exonaddresscts2_OUTFILE)s | python addressct2.py %(exonaddressccts2_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

	# apply H=1 filter to eliminate UMIs that arise from sequencing error for whole genes
	exoncfilt_OUTFILE = outdir+'/'+prefix+'.exon_address.cfilt.txt'
	cmd = 'cat %(exonaddressccts2_OUTFILE)s | python callfilter.py %(exoncfilt_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

	# clean up temporary files
	cmd = 'rm %(exonaddresscts2_OUTFILE)s %(exonaddressccts2_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

	# compute summary statistics on UMI-filtered data for exons (this time after cell barcode collapse)
	exoncfiltpdf_OUTFILE = outdir+'/'+prefix+'.exon_address.cfilt.pdf'
	exoncfiltchist_OUTFILE = outdir+'/'+prefix+'.exon_address.cfilt.chist.txt'
	cmd = 'python dsstats.py %(exoncfilt_OUTFILE)s %(exoncfiltpdf_OUTFILE)s %(exoncfiltchist_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

elif technology == 'PearSeq' or technology in ['CiteSeq5v2','CiteSeq3v3','CiteSeqTSB']:
	if ui.post_swap:
		pearaddress_INFILE=outdir+'/'+prefix+'.address.swap.decode.txt'
	tmp_OUTFILE = outdir+'/'+prefix+'.tmp'
	pearaddresscts_OUTFILE = outdir+'/'+prefix+'.pear_address.cts.txt'
	cmd = 'zcat %(pearaddress_INFILE)s.gz | python addressct.py %(tmp_OUTFILE)s %(pearaddresscts_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)
	
	pearaddressccts_OUTFILE = outdir+'/'+prefix+'.pear_address.ccts.txt'
	pear_trackbcs1_OUTFILE = outdir+'/'+prefix+'.pear_address.trackbcs1.txt'
	cmd = 'python collapse.py %(pearaddresscts_OUTFILE)s %(pearaddressccts_OUTFILE)s %(pear_trackbcs1_OUTFILE)s %(technology)s' % vars()
	print(cmd)
	os.system(cmd)	

	pearaddresscccts_OUTFILE = outdir+'/'+prefix+'.pear_address.fcts.txt'
	cmd = 'cat %(pearaddressccts_OUTFILE)s | python addressct2.py %(pearaddresscccts_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)
	
	cmd = 'rm %(pearaddressccts_OUTFILE)s %(pear_trackbcs1_OUTFILE)s %(pearaddresscts_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)

	pearcfiltpdf_OUTFILE = outdir+'/'+prefix+'.pear_address.fcts.pdf'
	pearcfiltchist_OUTFILE = outdir+'/'+prefix+'.pear_address.fcts.chist.txt'
	cmd = 'python dsstats.py %(pearaddresscccts_OUTFILE)s %(pearcfiltpdf_OUTFILE)s %(pearcfiltchist_OUTFILE)s' % vars()
	print(cmd)
	os.system(cmd)






