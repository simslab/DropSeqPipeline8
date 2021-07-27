# DropSeqPipeline8

This is a data processing pipeline for analyzing microwell- or DropSeq-like scRNA-seq data.  It can also be used for analyzing pooled plate-based scRNA-seq (as in Snyder et al, Science Immunology, 2019, for example). The pipeline in is written in Python, uses STAR for alignment, features automatic backup to an Amazon Web Services S3 bucket, and has the following dependencies:

-Python >=3.6
-bx.intervals package (https://bx-python.readthedocs.io/en/latest/lib/bx.intervals.html)
-Numpy 
-matplotlib
-STAR >=2.7.0d
-samtools >=1.9

Basic usage for analyzing microwell- or DropSeq-based scRNA-seq data is:

python DropSeqPipeline8.py -o output_directory -p output_prefix -r1 input_R1.fastq.gz -r2 input_R2.fastq.gz -r reference_genome_directory -g reference_genome_directory/annotation.gtf -v 65 -t 8 --s3path s3_bucket_name --technology DropSeqv2

where 65 is the sjdbOverhang parameter used to build the suffix array for STAR and 8 is the desired number of threads for parallelizing STAR and samtools.

