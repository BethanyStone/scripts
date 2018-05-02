#!/bin/bash
set -e
set -u

##### DNA METHYLATION RE-EXTRACTION #####

# script to re-extact DNA methylation from sam file output of bismark alignmnet (incase you've accidently deleted files you may be needing...)
# run in 4_bismark output sub-directory of wgbs_v1.0.sh workflow
# performs bismark_mtehylation_extractor and bismark2bedGraph commands to re-generate cytosine reports and context-specific bedGraph files...

### preparation ###

#required arguments
if [ "$#" -ne 1 ]; then
echo "USAGE: met-sign.sh <context> <file> <file path to met bedfile> <annotation file> <sample> <outname>"
exit 1
fi

# gather input variables
samfile=$1


### bismark methylation extraction ###

echo "Re-extracting cytosine report from SAM file..."

bismark_methylation_extractor --comprehensive --cytosine_report --CX --genome_folder ~/TAIR10_bs/  --report --buffer_size 8G -s ${samfile}


### bedgraph re-creation ###

echo "Done methylation extraction... Re-creating bedgraphs..."

bismark2bedGraph --CX CpG* -o ${samfile::-4}_CpG.bed
bismark2bedGraph --CX CHG* -o ${samfile::-4}_CHG.bed
bismark2bedGraph --CX CHH* -o ${samfile::-4}_CHH.bed

echo "DONE"
