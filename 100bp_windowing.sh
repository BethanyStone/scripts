#!/bin/bash

##### 100BP WINDOW GENERATION #####

### Requirements
# bedtools v2.25.0
# awk

# this script takes .bed.gz.bismark.cov.gz file (coverage file output from bismark alignment) and generates 100bp sliding windows for which methylation levels are summed and averaged across
# need to use windowing_template.sh to generate a "template" file (in bed format) containing windows of interest based on chromosome sizes (or whatever region you're interested in)
# run in directory containing .bed.gz.bismark.cov.gz file
# .bed.gz.bismark.cov.gz file contains columns as follows:
#	chromosome	Start	Stop	%methylated	count_unmethylated	count_methylated
# windowing template file contains columsn as follows:
#	chromosome	window_start	window_stop
# script adapted from https://www.biostars.org/p/70577/

### usage
if [ "$#" -ne 2 ]; then
echo "USAGE: 100bp_windowing.sh <bed file> <path to template windowing file"
echo "EXAMPLE: 100bp_windowing.sh Exp532_Col-0_1_r1_CHH.bed.gz.bismark.cov.gz ~/annotation/TAIR10/TAIR10_100bp_windows.bed"
exit 1
fi


### define variables
bed_file=$1		# bismark bed.cov file to be windowed
window_file=$2		# path to windowing template file 


### prepare coverage bed file

echo "Pre-sorting coverage bed file..."
# sortBed to sort coverage bed file
# sorts by chromosome and start position prior to bedtools intersect command to mmake process less memory intensive
#	-i = input bed file
sortBed -i ${bed_file} > ${bed_file}.sorted


### perform windowing

echo "Performing windowing..."

# intersect bed file with windowing file using bedtools intersect
# screens for overlaps between template windowing file and sorted coverage bed file
# each feature in A is compared to B in search of overlaps (order is important!)
#	-wa = write original entry in A for each overlap
#	-wb = write original entry in B for each overlap
#	-a = BAM/BED/GFF/VCF file
#	-b = BAM/BED/GFF/VCF file for which overlaps with A will be searched for
bedtools intersect -wa -wb -a ${window_file} -b ${bed_file}.sorted > ${bed_file}.int

# remove sorted coverage bed file
rm ${bed_file}.sorted


### summarise methylation metrics across windows

echo "Summarising methylation across windows..."

# sortBed to sort intersect bed file
#sort intersected bed file by chromosome and start position prior to input into bedtools merge command
sortBed -i ${bed_file}.int > ${bed_file}.int.sorted

# remove intersect bed file
rm ${bed_file}.int

# subset to columns of interest
# use awk to remove unneeded columns in sorted intersect bed file - bedtools intersect just pastes all columns in A against all columns in B for which there are overlaps
# want to keep columns as follows:
#	Chromosome	Start	Stop	Count_methylated	Count_unmethylated
awk -F '\t' '{print $1,$2,$3,$8,$9}' OFS='\t' ${bed_file}.int.sorted > ${bed_file}.int.sorted.sub

# remove sorted intersect bed file
rm ${bed_file}.int.sorted

# summarise methylation across windows
# bedtools groupby computes summary statistics on a column (-c) based on appropriate grouping column (-g)
# need bed file to be sorted, as when a change in grouping columns are detected a new summary is started
#	-i = input bed file
#	-g = column that will be used to group the input
#	-o = action to be applied to columns being summarised
bedtools groupby -i ${bed_file}.int.sorted.sub -g 1,2,3 -c 4,5 -o sum > ${bed_file}.grouped

# cant perform multiple actions in one groupby call, so need to compute number of methylation sites (just count of lines in bed file) separately...
bedtools groupby -i ${bed_file}.int.sorted.sub -g 1,2,3 -c 4 -o count > ${bed_file}.grouped2

bedtools map -a ${bed_file}.grouped -b ${bed_file}.grouped2 -c 4 > ${bed_file}.grouped.out

# remove subsetted sorted intesect bed file and intermediate grouped files
rm ${bed_file}.int.sorted.sub
rm ${bed_file}.grouped
rm ${bed_file}.grouped2

### format grouped.out bed file

echo "Tidying window file..."

# want final windowed bed file with following columns:
#	chromosome
#	window start
#	window stop
#	count of unmethylated sites across window (from reads underlying window)
#	count of methylated sites across window (from reads underlying window)
#	weighted average methylation of window (calculated as % of all methylation sites that are methylated in underlying reads)
#	number of methylation sites in window

# gather info for columns as above
#	-F = field separator ('\t' = tab-delimited)
#	OFS = output field separator ('\t')
#	$7 = percentage of methylated sites ((methylated/methylated + unmethylated)* 100)
awk -F '\t' '{$7=($4/($4+$5))*100;} {print $1,$2,$3,$4,$5,$7,$6}' OFS='\t' ${bed_file}.grouped.out > ${bed_file::-22}_100bp.bed

# remove grouped.out bed file
rm ${bed_file}.grouped.out

# gzip 100pb window output file
if [[ ${bed_file::-22}_100bp.bed != *gz ]];then
	gzip ${bed_file::-22}_100bp.bed;
fi

