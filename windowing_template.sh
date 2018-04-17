##### windowing_template.sh #####

#!/bin/bash

# commands used to generate a windowing template file for use in 100bp_windowing.sh
# takes a reference genome of interest and extracts chromosome sizes (or whatever other structures of interest), and uses this to create a bed file of sliding (non-overlapping) windows 
# script not set up to run automatically - need to manually run commands!




##### to generate GFF of chromosome sizes 
# this command is used to list name of chromosome (or scaffolds etc) and their size in basepairs
# must be a tab-delimited file with following columns:
#	chromosome	size
# based on https://www.biostars.org/p/173963/


### concatenate individual chromosome fasta files into single fasta file 
# only do this if needed!
# change fasta file names as needed
cat TAIR10_chr1.fa TAIR10_chr2.fa TAIR10_chr3.fa TAIR10_chr4.fa TAIR10_chr5.fa TAIR10_chrC.fa TAIR10_chrM.fa


### index reference genome fasta file
# samtools faidx command indexes a reference sequence in fasta format
# input single reference genome fasta file 
# change fasta file name as needed
samtools faidx TAIR10.fa


### pull out relevent information from indexed fasta file
cut -f1,2 TAIR10.fa.fai > TAIR10.chrom.sizes



##### generate windowing template file
### 100bp sliding windows
# takes chromosome size file and splits this into 100bp windows (or whatever size is specified)
# bedtools makewindows command is used:
#	-g = file containing chromosome sizes (in GFF)
#	-w = window size in basepairs
#	-s = specifies size of "slide" (i.e. where next window starts compared to first position in previous window) in basepairs
# change file names as needed
# ouput is a bed file with columns as follows:
#	chromsome	window_start	window_stop
bedtools makewindows -g TAIR10.chrom.sizes -w 99 -s 100 > TAIR10_100bp_windows.bed

