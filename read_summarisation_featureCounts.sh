##### read_summarisation_featureCounts.sh #####     

#!/bin/bash

set -e
set -u 


# this script performs read summarisation over features of interest following alignment to a refe$
# featureCounts is used to perform this (part of subread package)
# featureCounts documentation at: http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf
# required input files are as follows:
#       aligned reads in SAM or BAM format
#       a list of genomic features of interest in GTF of SAF format
#               chromosome names in annotation file MUST match chromosome names in reference geno$
#               SAF format annotation as follows: GeneID Chr Start End Strand
# SAF format annotation files can be created from bed file annotation using make_SAF_annotation.s$
# based on pedrocrisp's feeatureCounts summarisation and wgbs_pipeline_v1.0.sh


### usage ### 

if [ "$#" -lt 6 ]; then
echo "Missing required arguments!"
echo "USAGE: read_summarisation_featureCounts.sh <input bam file> <-se or -pe> <annotation_file.b$
echo "EXAMPLE: Read_summarisation_featureCounts.sh Exp532_Col-0_2_r1_sorted.bam -se /home/bethany$
exit 1
fi


### requirements ### 
bam=$1  # input bam file to be summarised
type=$2 # specifies whether aligned reads are single-end or paired-end reads
anno=$3 # path to annotation file
format=$4       # file format of annotation file provided (-gtf or -saf
fileID=$5       # fileID for output files
feat=$6         # feature reads are being summariesed across


### if SAF format annotation file is provided ###
# uses SAF formatted file for alignments

if [[ ${format} == "-saf" ]]; then

echo "Using SAF format annotation file..."

### read summarisation using featureCounts ###
# read summarisation over genomic features is performed using featureCounts
#       assigns reads by comparing mapping location of each base in the read with the genomic reg$
#       takes account of gaps line insertions, deletions, exon-exon junctions, structural variant$
#       hit is called if any overlap is found between the read and the feature
# arguments are as follows (only some of them, see documentation for full details):
#       -a = name of annotation file
#       -C = excludes chimeric fragments (two ends align to different chromosomes) from being cou$
#       -d = minimum fragment/template length to be considered (50 by default)
#       -D = maximum fragment/template length to be considered (50 by default)
#       -f = summarise reads at feature level (exon); summarised at meta-feature level by default$
#       -F = format of the annotation file ('GTF' or 'SAF')
#       -g = atribute for grouping features into metafeatures when GTF annotation is provided (ge$
#       -M = if specified, multimapping reads will be counted
#       -O = allow reads overlapping multiple features to be counted
#       -p = fragments are counted instead of reads (applicable for paired-end reads)
#       -s = strand specific reads counting will be performed (0=unstranded, 1=stranded, 2=reverse)
#       -t = specify feature type for GTF files (exon by default)
#       -T = number of threads to be used

# for single-end reads
if [[ ${type} == "-se" ]]; then 
echo "Performing counts for single-end reads..."
featureCounts -a ${anno} -F 'SAF' -C -T 2 -s 0 -o "${fileID}_${feat}.counts" ${bam}
fi

# for paired-end reads
if [[ ${type} == "-pe" ]]; then 
echo "Performing counts for paired-end reads..."
featureCounts -a ${anno} -F 'SAF' -p -C -T 2 -s 0 -o "${fileID}_${feat}.counts" ${bam}
fi

fi


### if GTF format annotation file is provided ###
# uses GTF formatted file for alignments

if [[ ${format} == "-gtf" ]]; then 

echo "Using GTF format annotation file..."

### read summarisation using featureCounts ###

# for single-end reads
if [[ ${type} == "-se" ]]; then 
echo "Performing counts for single-end reads..."
featureCounts -a ${anno} -F 'GTF' -g gene_id -C -T 2 -s 0 -o "${fileID}_${feat}.counts" ${bam}
fi

# for paired-end reads
if [[ ${type} == "-pe" ]]; then 
echo "Performing counts for paired-end reads..."
featureCounts -a ${anno} -F 'GTF' -g gene_id -p -C -T 2 -s 0 -o "${fileID}_${feat}.counts" ${bam}
fi

fi




