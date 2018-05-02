##### wgs_pipeline_v1.0.sh #####

#!/bin/bash

set -e
set -u

# #!/bin/bash specifies that script is written in bash code
# set -e = exit pipeline immediately if one of the commands executed returns an error
# set -u = treat unset variables and parameters as an error when performing parameter expansion

# this script takes single or paired-end fastq file from whole-genome sequencing, performs quality checking steps (fastqc), quality trimming (trimgalore!), aligns readas to reterence genome and indexes the resulting output
# script requires a subread indexed reference genome (generated as follows)
# to index genome:
#	subread-buildindex -o TAIR10_subread_index TAIR10_Chr1.fasta TAIR10_Chr2.fasta TAIR10_Chr3.fasta TAIR10_Chr4.fasta TAIR10_Chr5.fasta TAIR10_ChrC.fasta TAIR10_ChrM.fasta
# script is based on wgbs_pipeline_v1.0.sh, RNA-Seq_pipeline_v1.0.sh, and dtrain gdna_wgs_pipe.v1.sh


### usage ###

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: wgs_pipeline_v1.0.sh <-se or -pe> <R1> <R2> <path to subread indexed genome> <fileID output>"
exit 1
fi


##############
# SINGLE END #
##############

### preparation ###

# confirm single-end alignment
if [ "$1" == "-se" ]; then

# required arguments
if [ "$#" -ne 4 ]; then
echo "Required arguments missing for single-end alignment!"
echo "USAGE: wgs_pipeline_v1.0.sh <-se> <in fastq R1> <path to subread indexed genome> <fileID for output files>"
echo "EXAMPLE: wgs_pipeline_v1.sh -se fastq_file_1.fastq.gz /home/bethanys/TAIR10/TAIR10_subread_index fastq_file_1"
exit 1
fi

# gather input variables
type=$1 # specifies single-end alignment (-se or -pe)
fq=$2; #input fastq file to be aligned
index=$3; #path to subread indexed reference genome
fileID=$4; # fileID for output files
dow=$(date +"%F-%H-%m-%S") # date/time of alignment (don't need to provide this - it is done automatically)

echo "##################"
echo "Performing single-end WG-Seq alignment with the following parameters:"
echo "Type: ${type}"
echo "Input Files: ${fq}"
echo "genome index: ${index}"
echo "Output ID: ${fileID}"
echo "Time of analysis: ${dow}"
echo "##################"

# make sample directory (will be make in directory where alignment is being run from)
mkdir ${fileID}_wgs_${dow}

# move input fastq file into sample directory
mv $fq ${fileID}_wgs_${dow}

# change directories to sample directory
cd ${fileID}_wgs_${dow}



### initial fastqc ###
# fastqc performs quality checking on fastq files before they are aligned to a reference genome (see RNA-Seq_pipeline_v1.0.sh for details)

# gzip fastq file if it isn't already (should be though)
if [[ ${fq} != *.gz ]];then
gzip ${fq}
fq="${fq}.gz"
fi

echo "Initial FASTQC..."

# make directory for initial fastqc output files (1_fastqc)
mkdir 1_fastqc

# run fastqc on input fastq file using fastqc command
fastqc -t 2 ${fq} 2>&1 | tee -a ${fileID}_logs_${dow}.log

# move fastqc output to 1_fastqc directory
mv ${fq%%.fastq*}_fastqc* 1_fastqc



### trimgalore ###

echo "Done... trimming... "

# make directory for trimgalore output (2_trimgalore)
mkdir 2_trimgalore

# move to 2_trimgalore directory
cd 2_trimgalore/

# run trimgalore command to perform trimming of input fastq file 
# uses default settings (see wgbs_pipeline_v1.0.sh) and dont gzip option
trim_galore --dont_gzip ../${fq} 2>&1 | tee -a ../${fileID}_logs_${dow}.log

# move back to sample directory
cd ../



### trimmed fastqc ###

echo "Done... trimmed FASTQC..."

# make directory for trimmed fastq output files (3_trimmed_fastqc)
mkdir 3_trimmed_fastqc

# run fastqc command on "_trimmed.fq" output from trimgalore command
fastqc -t 2 2_trimgalore/${fq%%.fastq*}_trimmed.fq 2>&1 | tee -a ${fileID}_logs_${dow}.log

# move trimmed fastqc output files to 3_trimmed_fastqc directory
mv 2_trimgalore/${fq%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc



### organise files & directories ###

echo "Done... cleaning..."

# make directory for input fastq file ("0_rawfastq")
mkdir 0_rawfastq

# move input fastq file to 0_rawfastq directory
mv $fq 0_rawfastq



### subread alignment ###

echo "Perfoming subread alignment..."

# make directory for subread align output files (4_subread-align)
mkdir 4_subread-align

# move trimmed_fq file to 4_subread-align directory for alignment
mv 2_trimgalore/${fq%%.fastq*}_trimmed.fq -t 4_subread-align/

# change directories to 4_subread-align
cd 4_subread-align/

# run subread alignment
# see RNA-Seq_pipeline_v1.0.sh for subread-align details
#	-T 4 = threads/CPUs used for mapping
#	-t 1 = type of input sequencing data (0 for RNA seq, 1 for genomic DNA-seq)
#	-i = path to indexed reference genome
#	-r = fastq input file
#	-o = fileID for output file
subread-align -T 4 -M 0 -t 1 -i ${index} -r ${fq%%.fastq*}_trimmed.fq -o "${fileID}.bam" 2>&1 | tee -a ../${fileID}_logs_${dow}.log

# gzip trimmed_fq file following alignment
if [[ $fq%%.fastq}* != *".gz" ]]; then gzip ${fq%%.fastq*}_trimmed.fq; fi


### clean files and directories ###

echo "Alignment compete... making sorted bam file with index"

# set name for temporary bam files
tmpbam="${fileID}.bam"

# set name for output bam files
outbam="${fileID}.sorted.bam"

# samtools sort to sort bam file by chromosomal position, and save the sorted file
samtools sort -m 2G ${tmpbam} -o $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

# samtools index to make index of sorted bam file
samtools index $outbam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

# remove temporary bam files
rm -v ${tmpbam}

# move trimmed.fq file back to 2_trimgalore directory following alignment
mv *trimmed.fq.gz ../2_trimgalore/

echo "Alignment complete!"

fi




####
# PAIRED END
####

if [ "$1" == "PE" ]; then

if [ "$#" -ne 5 ]; then
echo "Missing required arguments for paired-end!"
echo "USAGE: gdna_wgs_pipe.v1.sh <PE> <R1> <R2> <subread indexed genome> <fileID output>"
echo "EXAMPLE: gdna_wgs_pipe.v1.sh PE sample_R1.fastq sample_R2.fastq $HOME/TAIR10/chromosomes/TAIR10_subread_index sample-r1"
exit 1
fi

#gather input variables
type=$1
fq1=$2;
fq2=$3;
index=$4; #path to genome index
fileID=$5;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing paired-end RNA-seq alignment with the following parameters:"
echo "Type: $type"
echo "Input Files: $fq1 $fq2"
echo "genome index: $index"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo "##################"

# make sample work directory
mkdir ${fileID}_gDNA_wgs_${dow}
mv $fq1 ${fileID}_gDNA_wgs_${dow}
mv $fq2 ${fileID}_gDNA_wgs_${dow}
cd ${fileID}_gDNA_wgs_${dow}

if [[ $fq1 != *.gz ]];then
gzip $fq1
fq1="${fq1}.gz"
fi

if [[ $fq2 != *.gz ]];then
gzip $fq2
fq2="${fq2}.gz"
fi

# initial fastqc
mkdir 1_fastqc
fastqc -t 2 $fq1 $fq2 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq1%%.fastq*}_fastqc* 1_fastqc
mv ${fq2%%.fastq*}_fastqc* 1_fastqc

echo "Performing quality-based read trimming... "

mkdir 2_trimgalore
cd 2_trimgalore
trim_galore --dont_gzip --paired ../$fq1 ../$fq2 2>&1 | tee -a ${fileID}_logs_${dow}.log
cd ../

# repeat fastqc
echo "FASTQC ..."

# fastqc again
mkdir 3_trimmed_fastqc
fastqc -t 2 2_trimgalore/${fq1%%.fastq*}_trimmed.fq 2_trimgalore/${fq2%%.fastq*}_trimmed.fq 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 2_trimgalore/${fq1%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc
mv 2_trimgalore/${fq2%%.fastq*}_trimmed_fastqc* -t 3_trimmed_fastqc

mkdir 0_fastq
mv $fq1 0_fastq
mv $fq2 0_fastq

# subread align
mkdir 4_subread-align
mv 2_trimgalore/${fq1%%.fastq*}_trimmed.fq -t 4_subread-align/
mv 2_trimgalore/${fq2%%.fastq*}_trimmed.fq -t 4_subread-align/
cd 4_subread-align/

echo "Beginning alignment ..."

# -t 0 = RNA-seq -t 1 = genomic DNA seq

subread-align -T 4 -t 1 -i ${index} -r ${fq1%%.fastq*}_trimmed.fq -R ${fq2%%.fastq*}_trimmed.fq -o "${fileID}.bam" 2>&1 | tee -a ../${fileID}_logs_${dow}.log

echo "cleaning..."

if [[ $fq1%%.fastq}* != *".gz" ]]; then gzip ${fq1%%.fastq*}_trimmed.fq; fi
if [[ $fq2%%.fastq}* != *".gz" ]]; then gzip ${fq2%%.fastq*}_trimmed.fq; fi

mv *trimmed.fq.gz ../2_trimgalore/

echo "Alignment complete"

fi

