#!/bin/bash

set -e 
set -u


##### CHECKING AVERAGE GENOME COVERAGE #####

# script adapted from dtrain/average_cov.sh
# computes average genome coverage of BS-seq data
# input file is sorted bismark bam output file - need to run script in directory containing this file...
# samtools depth command calculates read depth at each position (or region) of genome
# 	this is then averaged across the entire genome

samtools depth  *fq_bismark.sorted.sam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'

