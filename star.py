#! /usr/bin/python
import os
import sys

# The purpose of this program is to align single cell RNA-Seq data to the genome and transcriptome (Gencode annotations) using two-pass STAR.  

outdir = sys.argv[1]  # name of directory containing fastq (e.g. PJ018)
prefix = sys.argv[2]
fastq2_INFILE = sys.argv[3]  # name of fastq file
ref = sys.argv[4]
gtf = sys.argv[5]
threads = int(sys.argv[6])
overhang = int(sys.argv[7])
bam_OUTFILE = outdir+'/'+prefix+'.'

# assumes a maximum read length of 66 cycles
cmd = 'STAR --readFilesCommand zcat --genomeDir %(ref)s --sjdbOverhang %(overhang)d --sjdbGTFfile %(gtf)s --twopassMode Basic --runThreadN %(threads)d --readFilesIn %(fastq2_INFILE)s --outFileNamePrefix %(bam_OUTFILE)s --outSAMtype BAM Unsorted' % vars()
os.system(cmd)

