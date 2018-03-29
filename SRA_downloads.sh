############################# 
# SRA DOWNLOADS
#############################

# Script for performing bulk downloads of SRR files from SRA based on input list (as csv file) of SRRs
# Input csv file formatted as [Name, Rep, Run, Strategy, Instrument]
# Need to specify single end or paired end
# Run script from directory 
# Need to input either single end OR paired end reads - script not yet written to handle both at once

#!/bin/bash
set -eu

### usage
if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: SRA_download.sh <-se or -pe> <in csv_file>"
fi
	
	
### define arguments

strategy=$1		# paired end or single end reads?
csv=$2		# the input csv file of SRRs to be downloaded
dow=$(date +"%F-%H-%m-%S")		# timestamp


### make SRR directory

mkdir SRR_downloads_${dow}	# make directory to store downloads
cd SRR_downloads_${dow}		# move to directory to store downlads



##### SINGLE END #####
# Download SRRs, rename based on input csv file, and perform fastq-dump
# SRR files will be deleted, fastq files will be zipped
# SRR download based on col 5 of input csv file
# renames SRR file as Name.Rep.Strategy.Run so multiple runs can be concatenated together in later steps


### confirm single-end

if [ "$1" == "-se" ];then


### Download SRRs, rename, and perform fastq dump

tail -n +2 ../SRA_list.csv | tr -d '\r' | awk -F "," 'BEGIN {OFS="."} {print $1,$2,$4,$3,$5}' | while read i;		
do
	wget -nc -O ${i%.*} http://sra-download.ncbi.nlm.nih.gov/srapub/${i##*.};
	fastq-dump ${i%.*};
	rm ${i%.*};
done


### Concatenate multiple runs into single fastq file

runs=$(tail -n +2 ${csv} | tr -d '\r' | awk -F "," 'BEGIN {OFS="."} {print $1,$2,$4}' | uniq)

for i in $runs
do
	find ${i}* -exec cat {} > ${i}.fastq;
	gzip ${i}.fastq;
done


fi



##### PAIRED END #####


### Confirm paired ended

if [ "$1" == "-pe" ];then


### Download SRRs, rename, and perform fastq-dump
# Use --split-3 for paired end data

tail -n +2 ../SRA_list.csv | tr -d '\r' | awk -F "," 'BEGIN {OFS="."} {print $1,$2,$4,$3,$5}' | while read i;		
do
	wget -nc -O ${i%.*} http://sra-download.ncbi.nlm.nih.gov/srapub/${i##*.};
	fastq-dump --split-3 ${i%.*};
	rm ${i%.*};
	gzip ${i%.*}.fastq
done

# No commands as yet to concatenate paired end data containing multiple runs... (haven't come across any data that requires this yet)

fi
