##### make_SAF_annotation.sh #####
  
#!/bin/bash

set -e
set -u

# this script converts a bed file annotation into a SAF formatted annotation, which is required for edgeR DEG calling
# SAF file format contains following column
#       GeneID  Chr     Start   End     Strand
# input bed file annotation of interest into script to create the annotation
# NOTE: chromosome names in annotation file MUST match chromosome names in reference genome to which reads were aligned...
# run script in directory containing BED file annotation of interest

### usage ### 

if [ "$#" -lt 2 ]; then
echo "Missing required arguments!"
echo "USAGE: make_SAF_annotation.sh <input bam file> <fileID>"
echo "EXAMPLE: make_SAF_annotation.sh Araport11_gene.bed Araport11_gene"
exit 1
fi


### requirements ###
bed=$1  # bed file to be converted in SAF file
fileID=$2       # fileID for output SAF file


### make SAF file ###
awk -F'\t' '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' ${bed} > temp.saf
awk 'BEGIN {print "GeneID""\t""Chr""\t""Start""\t""End""\t""Strand"}{print}' temp.saf > temp2.saf

# remove temporary saf file
rm temp.saf

# rename chromosome names to match reference genome naming
 awk -F'\t' -vOFS='\t' '{ gsub("1", "Chr1", $2) ; gsub("2", "Chr2", $2) ; gsub("3", "Chr3", $2); gsub("4", "Chr4", $2); gsub("5", "Chr5", $2); gsub("^C$", "chloroplast", $2); gsub("^M$", "mitochondria", $2); print}' temp2.saf > ${fileID}.saf

# remove second temporary saf file
rm temp2.saf
