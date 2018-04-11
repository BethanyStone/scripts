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
# bedtools merge combines overlapping features (i.e windows) into a single feature spanning all of the combined features
# since windows are non-overlapping, methylation will be summarised for each window
#	-i = input file
#	-c = columns to be summarised
#	-o = action to be applied to columns being summarised
bedtools merge -i ${bed_file}.int.sorted.sub -c 4,5 -o sum,sum > ${bed_file}.merge

# remove subsetted sorted intesect bed file
rm ${bed_file}.int.sorted.sub


### format merged bed file

echo "Tidying window file..."

# want final windowed bed file with following columns:
#	chromosome
#	window start
#	window stop
#	count of unmethylated sites across window (from reads underlying window)
#	count of methylated sites across window (from reads underlying window)
#	average methylation of window (calculated as $ of all methylation sites that are methylated in underlying reads)
#	ideally, number of methylation sites in window (but haven't worked out how to calculate this properly yet...)

# gather info for columns as above
#	-F = field separator ('\t' = tab-delimited)
#	OFS = output field separator ('\t')
#	$6 = total # methylation sites in reads underlying window (count_methylated + count_unmethylated)
#	$7 = percentage of methylated sites (methylated/total sites * 100)
awk -F '\t' '{$6=$4+$5;} {$7=$4/$6*100;} {print $1,$2,$3,$4,$5,$7,$6}' OFS='\t' ${bed_file}.merge > ${bed_file}_100bp.bed

# remove merged bed file
rm ${bed_file}.merge

# gzip 100pb window output file
if [[ ${bed_file}_100bp.bed != *gz ]];then
	gzip ${bed_file}_100bp.bed;
fi

