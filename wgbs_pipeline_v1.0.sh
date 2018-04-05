##### wgbs_pipeline_v1.0.sh #####

#!/bin/bash
set -e 
set -u

# #!/bin/bash specifies that script is written in bash code
# set -e = exit pipeline immediately if one of the commands executed returns a non-sero status (if the pipeline has an error)
# set -u = treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion (an error message will be written to the standard error and a non-interactive shell will exit)

# this script takes a single-or paired-end fastq file, runs it through quality checking steps (fastqc), quality trimming (trimgalore), and performs a bismark alignment to align BS-reads to a reference genome and call methylated cytosines
# script generates per-cytosine bed files and 100bp window wig files for CpG, CHG and CHH methylation levels
# script is based on SRE Bisulfite sequence analysis pipelien and dtrain/NGS-scripts/wgbs_pipelinev0.4.sh

# requires bowtie indexed reference genome
# to index genome:
#	 bowtie1 = bismark_genome_preparation /path/to/genome
#	 bowtie2 = bismark_genome_preparation --bowtie2 /path/to/genome
# bowtie2 seems to be better at aligning shorter reads, compared to bowtie1 (tend to get more aligned reads with bowtie2 than bowtie1)
# will use bowtie2 in this script


### usage ###

if [ "$#" -lt 4 ]; then
echo "Missing required arguments!"
echo "USAGE: /home/bethany.stone/NGS-scripts/wgbs_pipelinev1.0.sh <-pe, -se, -se_epi, or -pese> <in fastq R1> <in fastq R2 (if -pe)> <path to bismark indexed genome folder> <fileID for output files>"
exit 1
fi




##############
# SINGLE END #
##############

### preparation ###

# confirm single-end alignment 
if [ "$1" == "-se" ];then

# required arguments
if [ "$#" -ne 4 ]; then
echo "Required arguments missing for single-end alignment!"
echo "USAGE: wgbs_pipelinev1.0.sh <-se> <in fastq R1> <path to bismark genome folder> <fileID for output files>"
echo "EXAMPLE: wgbs_pipelinev1.0.sh -se fastq_file_1.fastq.qz ~/reference_genomes/TAIR10 fastq_file_1"
exit 1
fi

# gather input variables
type=$1; # specifies single-end (-se) or paired-end (-pe) reads?
fq=$2; # input fastq file to be aligned
bisgenome=$3; # path to bismark indexed reference genome (give path to folder containing bismark genome prep files
fileID=$4; # fileID for output files
dow=$(date +"%F-%H-%m-%S") # date/time of alignment (don't need to provide this - it is added automatically)

echo "##################"
echo "Performing Bismark single-end alignment with the following parameters:"
echo "Type: ${type}"
echo "Input File: ${fq}"
echo "Path to bismark genome folder: ${bisgenome}"
echo "Output ID: ${fileID}"
echo "Time of analysis: ${dow}"
echo "Full logfile of steps: ${fileID}_${dow}.log"
echo "##################"

# make sample directory (will be made in directory where alignment is being run from)
mkdir ${fileID}_${dow}

# move input fastq file into sample directory
mv ${fq} ${fileID}_${dow}

# change directories to sample directory
cd ${fileID}_${dow}


### initial fastqc ###
# fastqc performs quality checking on fastq files before they are aligned to a reference genome (see RNA-Seq_pipeline_v1.0 for details)

echo "Initial FASTQC..."

# make directory for initial fastqc output files ("1_fastqc")
mkdir 1_fastqc

# run fastqc on input fastq file using fastqc command
# tee command adds output to log file (see RNA-Seq_pipeline_v1.0 for details)
fastqc ${fq} 2>&1 | tee -a ${fileID}_${dow}.log

# move fastqc output files to 1_fastqc directory
mv ${fq%%.fastq*}_fastqc* 1_fastqc #


### trimgalore ###

echo "Done... trimming..."

# trimgalore used Cutadapt and Fastqc to perform trimming of adapter sequences and poor quality reads from input fastq file prior to alignment to a reference genome
# uses first 13bp of Illumina standard adapters by default
# accepts gzipped fastq or unzipped fastq files
# trimgalore options (main ones):
#	-v = print version information and exit
#	-q = trim low-quality end reads in addition to adapter removal (default phred score is 20)
#	--phred33 = cutadapt uses ASCII+33 quality scores (default)
#	--phred65 = cutadapt uses ASCII+64 quality scores
#	--fastqc = run fastqc in default mode on fastq once trimming is complete
#	-a = adapter sequence to be trimmed (default is to autodetect if its Illumina universal, Nextera transposase or Illumina sRNA adapters)
#	-s = overlap with adapter sequence required to trim a sequence (default is 1- very stringent)
#	-e = max error rate allowed (no. errors/length of match; default is 0.1)
#	--gzip = compress output file with gzip
#	--length = discard reads shorter than length due to quailty/adapter trimming (default is 20bp)
#	-o = output written to this directory instead of current directory

# make directory for trimgalore output ("2_trimgalore")
mkdir 2_trimgalore

# move to 2_trimgalore directory
cd 2_trimgalore

# run trim_galore command to perform trimming of input fastq file
trim_galore ../${fq} 2>&1 | tee -a ../${fileID}_${dow}.log

# move back to sample directory
cd ../


### trimmed fastqc ###

echo "Done... trimmed FASTQC..."

# make directory for trimmed fastqc output files ("3_trimmed_fastqc")
mkdir 3_trimmed_fastqc

# run fastqc command on "_trimmed.fq" output file from trimgalore command
fastqc 2_trimgalore/${fq%%.fastq*}_trimmed.fq* 2>&1 | tee -a ${fileID}_${dow}.log

# move trimmed fastqc output files to 3_trimmed_fastqc directory
mv 2_trimgalore/${fq%%.fastq*}_trimmed_fastqc* 3_trimmed_fastqc


### organise files & directories ###

echo "Done... cleaning..."

# make directory for input fastq file ("0_fastq")
mkdir 0_rawfastq

# move input fastq file to 0_rawfastq directory
mv ${fq} 0_rawfastq


### bismark alignment ###

echo "Performing bismark alignment..."

# bismark aligns BS-seq data to a reference genome and performs cytosine methylation calling at the same time
# BS-seq reads are mapped using Bowtie1 or Bowtie2 short read aligner
# reference genome needs to be indexed using bismark_genome_preparation command, specifying bowtie2 option (see info at top of script)
# Bismark options when using bowtie2 (main ones):
#	-N = number of mismatches allowed in a seed alignment during multiseed alignment (0 or 1 where 1 is more sensitive but slower; default is 0)
#	-L = length of seed substrings to align during multiseed alignment; smaller values make alignments slower but more sensitive (default is 20, --sensitive preset)#	--sam = output written in SAM format instead of default BAM format (bismark will attempt to find samtools, or if samtools cannot be found SAM output will be compressed with gzip giving a sam.gz output file)

# make directory for bismark alignment output
mkdir 4_bismark_alignment

# move to 4_bismark)_alignment directory
cd 4_bismark_alignment

# to run bismark using bowtie1 (not used, but option there if needed - just comment out bowtie2 command):
#bismark --bowtie1 --sam -n 2 -l 20 ../../$genome_path ../2_trimgalore/${fq_file%%.fastq*}_trimmed.fq* 2>&1 | tee -a ../${fileID}_logs_${dow}.log

# command to run bismark using bowtie2 (default) on "_trimmed.fq" output file from trimgalore command
bismark --bowtie2 --sam -N 0 -L 20 ${bisgenome} ../2_trimgalore/${fq%%.fastq*}_trimmed.fq* 2>&1 | tee -a ../${fileID}_${dow}.log


### sam to bam file conversion ###

echo "Alignment complete... converting SAM to BAM..."

# samtools view to convert sam file (output from bismark alignment) to bam file
# samtools view options (main ones):
#	-b = output in BAM format
#	-C = output in CRAM format
#	-1 = enables fast BAM compression
#	-u = output uncompressed BAM
#	-h = include header in output
#	-o FILE = output to FILE
#	-S = input in SAM format (ignored for compatibility with previous samtools versions; correct format should be automatically detected by samtools)

samtools view -b -S -h ${fq%%.fastq*}_trimmed*.sam > ${fq%%.fastq*}_trimmed.fq_bismark.bam

# samtools sort command to sort bam file by chromosomal position
# samtools sort options:
#	-l = set compression level, from 0 (uncompressed) to 9 (best)
#	-m = set max memory per thread; suffix K/M/G recodnised
#	-n = sort by read name
#	-o = write final output to FILE rather than standard output
#	-T = write temporary files to PREFIX.nnn.bam
#	-@ = threads
#	-O = output format

samtools sort ${fq%%.fastq*}_trimmed.fq_bismark.bam ${fq%%.fastq*}_trimmed.fq_bismark.sorted 2>&1 | tee -a ../${fileID}_${dow}.log

# samtools index command to further compress bam file by indexing
# samtools index options:
#	-b = Generate BAI-format index for BAM files (default)
#	-c = Generate CSI-format index for BAM files
#	-m = min interval size for SCI indices (?)

samtools index ${fq%%.fastq*}_trimmed.fq_bismark.sorted.bam 2>&1 | tee -a ../${fileID}_${dow}.log


### bismark methylation extraction ###

echo "Sorting done... performing methylation extraction..."

# bismark_methylation_extractor operates on bismark result files and extracts methylation call for every cytosine analysed
# outputs a tab-delimited .txt file with:
#	1. seq-ID (CpG, CHG or CHH)
#	2. methylation state
#	3. chromosome
#	4. start position
#	5. methylation call
# bismark methylation_extractor options (only used options are listed):
#	--comprehensive = merve all for possible strand-specific methylation info into context-dependend output files (CpG, CHG and CHH)
#	--report = prints short methylation summary and parameters used to run script (default is on) 
#	-s = input files are bismark result files generated from single-end reads
#	-p = input files are bismark result files generated from paired-end reads
#	--buffer_size = specify main memory sort buffer for sorting methylation information (specify percentage, or can use K/G/M suffixes)

# run bismark_methylation_extractor command on .sam output file from bismark alignment
bismark_methylation_extractor --comprehensive --report --buffer_size 8G -s ${fq%%.fastq*}_trimmed*.sam 2>&1 | tee -a ../${fileID}_${dow}.log

# methylation extraction with full cytosine report (not used - uncomment to use if needed) 
##  bismark_methylation_extractor --comprehensive --cytosine_report --CX --genome_folder ~/TAIR10_bs/  --report --buffer_size 8G -s *.sam


### bedgraph creation ### 

echo "Creating bedGraphs"

# bismark2bedGraph outputs a bedGraph format of bismark methylation extractor methylation report (chromosome	start position	end position	methylation %)
# will be sorted by chromosomal coordinates 
# also adds coverate file that contains two additional columns of information about the coverate of reads at a cytosine position - count methylation and count unmethylated
# by default, only cytosines in CpG context are considered
# bismark2bedGraph options:
#	-- CX = extends command to consider CHG and CHH cytosine contexts

bismark2bedGraph --CX CpG* -o ${fileID}_CpG.bed 2>&1 | tee -a ../${fileID}_${dow}.log
bismark2bedGraph --CX CHG* -o ${fileID}_CHG.bed 2>&1 | tee -a ../${fileID}_${dow}.log
bismark2bedGraph --CX CHH* -o ${fileID}_CHH.bed 2>&1 | tee -a ../${fileID}_${dow}.log

# move back to sample directory
cd ../

# make directory for output files ("5_output_files") 
mkdir 5_output_files

# move .bed files from bismark2bedGraph commands to 5_output_files directory
mv 4_bismark_alignment/*.bed* 5_output_files


### 100bp window creation ###

# code from SRE that creates 100bp windows, and determines their methylation levels (for each sequence context)
# C_context_window_SREedits.pl won't work with perl v5.22.1... works with perl v5.22.1... can generate windowing data in a similar way using bedtools if needed.
# so, comment or uncomment this bit of code as needed

#perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CpG* 100 0 ${fileID}_CpG 2>&1 | tee -a ${fileID}_${dow}.log
#mv 4_bismark_alignment/CpG*.wig 5_output_files/${fileID}_CpG_100bp.wig

#perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CHG* 100 0 ${fileID}_CHG 2>&1 | tee -a ${fileID}_${dow}.log
#mv 4_bismark_alignment/CHG*.wig 5_output_files/${fileID}_CHG_100bp.wig
#perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CHH* 100 0 ${fileID}_CHH 2>&1 | tee -a ${fileID}_${dow}.log
#mv 4_bismark_alignment/CHH*.wig 5_output_files/${fileID}_CHH_100bp.wig


### gather pipeline metrics ###

echo "Providing pipeline metrics to wgbs pipeline logfile..."

# get all relevant numbers for final log summary
bismark_version=$(bismark --version | grep "Bismark Version:" | cut -d":" -f2 | tr -d ' ')
samtools_version=$(samtools 3>&1 1>&2 2>&3 | grep "Version:" | cut -d' ' -f2 | tr -d ' ')

map_ef=$(grep 'Mapping efficiency:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
unique_aln=$(grep 'Number of alignments with a unique best hit from the different alignments:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
no_aln=$(grep 'Sequences with no alignments under any condition:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
multi_aln=$(grep 'Sequences did not map uniquely:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
cpg_per=$(grep 'C methylated in CpG context:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chg_per=$(grep 'C methylated in CHG context:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chh_per=$(grep 'C methylated in CHH context:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)

if [[ ${fq} == *gz* ]];then
	raw_reads=$(zcat 0_rawfastq/*.gz | wc -l)
	raw_reads=$(($raw_reads / 4 ));fi

if [[ ${fq} != *gz* ]];then
	raw_reads=$(wc -l < 0_rawfastq/$fq_file)
	raw_reads=$(($raw_reads / 4 ));fi

if [[ ${fq} == *gz* ]];then
	flt_reads=$(zcat 2_trimgalore/*.gz | wc -l)
	flt_reads=$(($flt_reads / 4));fi

if [[ ${fq} != *gz* ]];then
	flt_reads=$(wc -l < 2_trimgalore/${fq%%.fastq*}_trimmed.fq*)
	flt_reads=$(($flt_reads / 4));fi
fi


# add it to the full pipeline logfile
printf "${dow}\t${fq_file}\t${fileID}\t${genome_path##../}\t${type:1}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n"
printf "${dow}\t${fq_file}\t${fileID}\t${genome_path##../}\t${type:1}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n" >> $HOME/wgbs_se_pipeline_analysis_record.log


# compress sam and unsorted bam files

echo "Compressing files..."

find -name "*.sam" | xargs pigz
find -name "*bismark.bam" | xargs rm

fi



######################################
# SINGLE END EPIGENOME
######################################

All this does differently is trip the first 6bp from the fastq reads as per epicentre protocol to eliminate issues of random hexamer binding... (difference is in the TrimGalore! command)


### preparation ###

# confirm single-end epigenome alignment (for bisulfite-sequencing data)
if [ "$1" == "-se_epi" ];then

# required arguments
if [ "$#" -ne 4 ]; then
echo "Required arguments missing for single-end alignment!"
echo "USAGE: wgbs_pipelinev1.0.sh <-se_epi> <in fastq R1> <path to bismark genome folder> <fileID for output files>"
echo "EXAMPLE: wgbs_pipelinev1.0.sh -se_ep BS-seq_fastq_file_1.fastq.qz ~/reference_genomes/TAIR10 fastq_file_1"
exit 1
fi

# gather input variables
type=$1; # specifies single-end (-se) or paired-end (-pe) reads?
fq=$2; # input fastq file to be aligned
bisgenome=$3; # path to bismark indexed reference genome (give path to folder containing bismark genome prep files
fileID=$4; # fileID for output files
dow=$(date +"%F-%H-%m-%S") # date/time of alignment (don't need to provide this - it is added automatically)

echo "##################"
echo "Performing Bismark single-end alignment with the following parameters:"
echo "Type: ${type}"
echo "Input File: ${fq}"
echo "Path to bismark genome folder: ${bisgenome}"
echo "Output ID: ${fileID}"
echo "Time of analysis: ${dow}"
echo "Full logfile of steps: ${fileID}_${dow}.log"
echo "##################"

# make sample directory (will be made in directory where alignment is being run from)
mkdir ${fileID}_${dow}

# move input fastq file into sample directory
mv ${fq} ${fileID}_${dow}

# change directories to sample directory
cd ${fileID}_${dow}


### initial fastqc ###
# fastqc performs quality checking on fastq files before they are aligned to a reference genome (see RNA-Seq_pipeline_v1.0 for details)

echo "Initial FASTQC..."

# make directory for initial fastqc output files ("1_fastqc")
mkdir 1_fastqc

# run fastqc on input fastq file using fastqc command
# tee command adds output to log file (see RNA-Seq_pipeline_v1.0 for details)
fastqc ${fq} 2>&1 | tee -a ${fileID}_${dow}.log

# move fastqc output files to 1_fastqc directory
mv ${fq%%.fastq*}_fastqc* 1_fastqc #


### trimgalore ###

echo "Done... trimming..."

# trimgalore used Cutadapt and Fastqc to perform trimming of adapter sequences and poor quality reads from input fastq file prior to alignment to a reference genome
# uses first 13bp of Illumina standard adapters by default
# accepts gzipped fastq or unzipped fastq files
# trimgalore options (main ones):
#	-v = print version information and exit
#	-q = trim low-quality end reads in addition to adapter removal (default phred score is 20)
#	--phred33 = cutadapt uses ASCII+33 quality scores (default)
#	--phred65 = cutadapt uses ASCII+64 quality scores
#	--fastqc = run fastqc in default mode on fastq once trimming is complete
#	-a = adapter sequence to be trimmed (default is to autodetect if its Illumina universal, Nextera transposase or Illumina sRNA adapters)
#	-s = overlap with adapter sequence required to trim a sequence (default is 1- very stringent)
#	-e = max error rate allowed (no. errors/length of match; default is 0.1)
#	--gzip = compress output file with gzip
#	--length = discard reads shorter than length due to quailty/adapter trimming (default is 20bp)
#	-o = output written to this directory instead of current directory
#	--clip_R1 6 = trim first 6 reads from read 1 (single ended reads anyway) to remove issues with random hexamer binding

# make directory for trimgalore output ("2_trimgalore")
mkdir 2_trimgalore

# move to 2_trimgalore directory
cd 2_trimgalore

# run trim_galore command to perform trimming of input fastq file
trim_galore --clip_R1 6 ../$fq_file 2>&1 | tee -a ../${fileID}_${dow}.log

# move back to sample directory
cd ../


### trimmed fastqc ###

echo "Done... trimmed FASTQC..."

# make directory for trimmed fastqc output files ("3_trimmed_fastqc")
mkdir 3_trimmed_fastqc

# run fastqc command on "_trimmed.fq" output file from trimgalore command
fastqc 2_trimgalore/${fq%%.fastq*}_trimmed.fq* 2>&1 | tee -a ${fileID}_${dow}.log

# move trimmed fastqc output files to 3_trimmed_fastqc directory
mv 2_trimgalore/${fq%%.fastq*}_trimmed_fastqc* 3_trimmed_fastqc


### organise files & directories ###

echo "Done... cleaning..."

# make directory for input fastq file ("0_fastq")
mkdir 0_rawfastq

# move input fastq file to 0_rawfastq directory
mv ${fq} 0_rawfastq


### bismark alignment ###

echo "Performing bismark alignment..."

# bismark aligns BS-seq data to a reference genome and performs cytosine methylation calling at the same time
# BS-seq reads are mapped using Bowtie1 or Bowtie2 short read aligner
# reference genome needs to be indexed using bismark_genome_preparation command, specifying bowtie2 option (see info at top of script)
# Bismark options when using bowtie2 (main ones):
#	-N = number of mismatches allowed in a seed alignment during multiseed alignment (0 or 1 where 1 is more sensitive but slower; default is 0)
#	-L = length of seed substrings to align during multiseed alignment; smaller values make alignments slower but more sensitive (default is 20, --sensitive preset)#	--sam = output written in SAM format instead of default BAM format (bismark will attempt to find samtools, or if samtools cannot be found SAM output will be compressed with gzip giving a sam.gz output file)

# make directory for bismark alignment output
mkdir 4_bismark_alignment

# move to 4_bismark)_alignment directory
cd 4_bismark_alignment

# to run bismark using bowtie1 (not used, but option there if needed - just comment out bowtie2 command):
#bismark --bowtie1 --sam -n 2 -l 20 ../../$genome_path ../2_trimgalore/${fq_file%%.fastq*}_trimmed.fq* 2>&1 | tee -a ../${fileID}_logs_${dow}.log

# command to run bismark using bowtie2 (default) on "_trimmed.fq" output file from trimgalore command
bismark --bowtie2 --sam -N 0 -L 20 ${bisgenome} ../2_trimgalore/${fq%%.fastq*}_trimmed.fq* 2>&1 | tee -a ../${fileID}_${dow}.log


### sam to bam file conversion ###

echo "Alignment complete... converting SAM to BAM..."

# samtools view to convert sam file (output from bismark alignment) to bam file
# samtools view options (main ones):
#	-b = output in BAM format
#	-C = output in CRAM format
#	-1 = enables fast BAM compression
#	-u = output uncompressed BAM
#	-h = include header in output
#	-o FILE = output to FILE
#	-S = input in SAM format (ignored for compatibility with previous samtools versions; correct format should be automatically detected by samtools)

samtools view -b -S -h ${fq%%.fastq*}_trimmed*.sam > ${fq%%.fastq*}_trimmed.fq_bismark.bam

# samtools sort command to sort bam file by chromosomal position
# samtools sort options:
#	-l = set compression level, from 0 (uncompressed) to 9 (best)
#	-m = set max memory per thread; suffix K/M/G recodnised
#	-n = sort by read name
#	-o = write final output to FILE rather than standard output
#	-T = write temporary files to PREFIX.nnn.bam
#	-@ = threads
#	-O = output format

samtools sort ${fq%%.fastq*}_trimmed.fq_bismark.bam ${fq%%.fastq*}_trimmed.fq_bismark.sorted 2>&1 | tee -a ../${fileID}_${dow}.log

# samtools index command to further compress bam file by indexing
# samtools index options:
#	-b = Generate BAI-format index for BAM files (default)
#	-c = Generate CSI-format index for BAM files
#	-m = min interval size for SCI indices (?)

samtools index ${fq%%.fastq*}_trimmed.fq_bismark.sorted.bam 2>&1 | tee -a ../${fileID}_${dow}.log


### bismark methylation extraction ###

echo "Sorting done... performing methylation extraction..."

# bismark_methylation_extractor operates on bismark result files and extracts methylation call for every cytosine analysed
# outputs a tab-delimited .txt file with:
#	1. seq-ID (CpG, CHG or CHH)
#	2. methylation state
#	3. chromosome
#	4. start position
#	5. methylation call
# bismark methylation_extractor options (only used options are listed):
#	--comprehensive = merve all for possible strand-specific methylation info into context-dependend output files (CpG, CHG and CHH)
#	--report = prints short methylation summary and parameters used to run script (default is on) 
#	-s = input files are bismark result files generated from single-end reads
#	-p = input files are bismark result files generated from paired-end reads
#	--buffer_size = specify main memory sort buffer for sorting methylation information (specify percentage, or can use K/G/M suffixes)

# run bismark_methylation_extractor command on .sam output file from bismark alignment
bismark_methylation_extractor --comprehensive --report --buffer_size 8G -s ${fq%%.fastq*}_trimmed*.sam 2>&1 | tee -a ../${fileID}_${dow}.log

# methylation extraction with full cytosine report (not used - uncomment to use if needed) 
##  bismark_methylation_extractor --comprehensive --cytosine_report --CX --genome_folder ~/TAIR10_bs/  --report --buffer_size 8G -s *.sam


### bedgraph creation ### 

echo "Creating bedGraphs"

# bismark2bedGraph outputs a bedGraph format of bismark methylation extractor methylation report (chromosome	start position	end position	methylation %)
# will be sorted by chromosomal coordinates 
# also adds coverate file that contains two additional columns of information about the coverate of reads at a cytosine position - count methylation and count unmethylated
# by default, only cytosines in CpG context are considered
# bismark2bedGraph options:
#	-- CX = extends command to consider CHG and CHH cytosine contexts

bismark2bedGraph --CX CpG* -o ${fileID}_CpG.bed 2>&1 | tee -a ../${fileID}_${dow}.log
bismark2bedGraph --CX CHG* -o ${fileID}_CHG.bed 2>&1 | tee -a ../${fileID}_${dow}.log
bismark2bedGraph --CX CHH* -o ${fileID}_CHH.bed 2>&1 | tee -a ../${fileID}_${dow}.log

# move back to sample directory
cd ../

# make directory for output files ("5_output_files") 
mkdir 5_output_files

# move .bed files from bismark2bedGraph commands to 5_output_files directory
mv 4_bismark_alignment/*.bed* 5_output_files


### 100bp window creation ###

# code from SRE that creates 100bp windows, and determines their methylation levels (for each sequence context)
# C_context_window_SREedits.pl won't work with perl v5.22.1... works with perl v5.22.1... can generate windowing data in a similar way using bedtools if needed.
# so, comment or uncomment this bit of code as needed

#perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CpG* 100 0 ${fileID}_CpG 2>&1 | tee -a ${fileID}_${dow}.log
#mv 4_bismark_alignment/CpG*.wig 5_output_files/${fileID}_CpG_100bp.wig

#perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CHG* 100 0 ${fileID}_CHG 2>&1 | tee -a ${fileID}_${dow}.log
#mv 4_bismark_alignment/CHG*.wig 5_output_files/${fileID}_CHG_100bp.wig
#perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CHH* 100 0 ${fileID}_CHH 2>&1 | tee -a ${fileID}_${dow}.log
#mv 4_bismark_alignment/CHH*.wig 5_output_files/${fileID}_CHH_100bp.wig


### gather pipeline metrics ###

echo "Providing pipeline metrics to wgbs pipeline logfile..."

# get all relevant numbers for final log summary
bismark_version=$(bismark --version | grep "Bismark Version:" | cut -d":" -f2 | tr -d ' ')
samtools_version=$(samtools 3>&1 1>&2 2>&3 | grep "Version:" | cut -d' ' -f2 | tr -d ' ')

map_ef=$(grep 'Mapping efficiency:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
unique_aln=$(grep 'Number of alignments with a unique best hit from the different alignments:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
no_aln=$(grep 'Sequences with no alignments under any condition:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
multi_aln=$(grep 'Sequences did not map uniquely:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
cpg_per=$(grep 'C methylated in CpG context:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chg_per=$(grep 'C methylated in CHG context:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chh_per=$(grep 'C methylated in CHH context:' 4_bismark_alignment/${fq%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)

if [[ ${fq} == *gz* ]];then
	raw_reads=$(zcat 0_rawfastq/*.gz | wc -l)
	raw_reads=$(($raw_reads / 4 ));fi

if [[ ${fq} != *gz* ]];then
	raw_reads=$(wc -l < 0_rawfastq/$fq_file)
	raw_reads=$(($raw_reads / 4 ));fi

if [[ ${fq} == *gz* ]];then
	flt_reads=$(zcat 2_trimgalore/*.gz | wc -l)
	flt_reads=$(($flt_reads / 4));fi

if [[ ${fq} != *gz* ]];then
	flt_reads=$(wc -l < 2_trimgalore/${fq%%.fastq*}_trimmed.fq*)
	flt_reads=$(($flt_reads / 4));fi
fi


# add it to the full pipeline logfile
printf "${dow}\t${fq_file}\t${fileID}\t${genome_path##../}\t${type:1}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n"
printf "${dow}\t${fq_file}\t${fileID}\t${genome_path##../}\t${type:1}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n" >> $HOME/wgbs_se_pipeline_analysis_record.log


# compress sam and unsorted bam files

echo "Compressing files..."

find -name "*.sam" | xargs pigz
find -name "*bismark.bam" | xargs rm

fi




######################################
# PAIRED END
######################################

if [ "$1" == "-pe" ];then

if [ "$#" -ne 5 ]; then
echo "Missing required arguments for paired-end!"
echo "USAGE: wgbs_se_pipelinev0.2.sh <-pe> <in fastq R1> <in fastq R2> <path to bismark genome folder> <fileID for output files>"
exit 1
fi
#gather input variables
type=$1; #identifying paired end or single end mode
fq_file1=$2; #the input R1 fastq reads
fq_file2=$3; #R2 reads
genome_path=$4; #the path to the genome to be used (bismark genome prepped)
fileID=$5;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing Bismark paired-end alignment with the following parameters:"
echo "Type: $type"
echo "Input File: $fq_file1"
echo "Path to bismark genome folder: $genome_path"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo ""
echo "Full logfile of steps: ${fileID}_logs_${dow}.log"
echo "##################"

#develop directory tree
mkdir ${fileID}_wgbspipeline_${dow}
mv $fq_file1 ${fileID}_wgbspipeline_${dow}
mv $fq_file2 ${fileID}_wgbspipeline_${dow}
cd ${fileID}_wgbspipeline_${dow}

#fastqc
mkdir 1_fastqc
fastqc $fq_file1 $fq_file2 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq_file1%%.fastq*}_fastqc* 1_fastqc #
mv ${fq_file2%%.fastq*}_fastqc* 1_fastqc #


#trim_galore
mkdir 2_trimgalore
cd 2_trimgalore
trim_galore --paired ../$fq_file1 ../$fq_file2 2>&1 | tee -a ../${fileID}_logs_${dow}.log
cd ../

#fastqc_again
mkdir 3_trimmed_fastqc
fastqc 2_trimgalore/${fq_file1%%.fastq*}_val_1.fq* 2_trimgalore/${fq_file2%%.fastq*}_val_2.fq* 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 2_trimgalore/${fq_file1%%.fastq*}_val_1_fastqc* 3_trimmed_fastqc
mv 2_trimgalore/${fq_file2%%.fastq*}_val_2_fastqc* 3_trimmed_fastqc


mkdir 0_rawfastq
mv $fq_file1 0_rawfastq
mv $fq_file2 0_rawfastq

#bismark
mkdir 4_bismark_alignment
cd 4_bismark_alignment
#bismark --bowtie1 --sam -n 2 -l 20 ../../$genome_path -1 ../2_trimgalore/${fq_file1%%.fastq*}_val_1.fq* -2 ../2_trimgalore/${fq_file2%%.fastq*}_val_2.fq* 2>&1 | tee -a ../${fileID}_logs_${dow}.log
bismark	--bowtie2 --sam -N 0 -L 20 ../../$genome_path -1 ../2_trimgalore/${fq_file1%%.fastq*}_val_1.fq* -2 ../2_trimgalore/${fq_file2%%.fastq*}_val_2.fq* 2>&1 | tee -a ../${fileID}_logs_${dow}.log

#sam to bam
samtools view -b -S -h ${fq_file1%%.fastq*}_val_1.fq*_bismark_pe.sam > ${fq_file1%%.fastq*}_val_1.fq_bismark_pe.bam
samtools sort ${fq_file1%%.fastq*}_val_1.fq_bismark_pe.bam ${fq_file1%%.fastq*}_val_1.fq_bismark_pe.sorted 2>&1 | tee -a ../${fileID}_logs_${dow}.log
samtools index ${fq_file1%%.fastq*}_val_1.fq_bismark_pe.sorted.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
#sam sort for MethylKit
#grep -v '^[[:space:]]*@' ${fq_file1%%.fastq*}_val_1.fq_bismark_pe.sam | sort -k3,3 -k4,4n > ${fq_file1%%.fastq*}_val_1.fq_bismark_pe.sorted.sam
#methylation extraction PAIRED END (-p)
bismark_methylation_extractor --comprehensive --report --buffer_size 8G -p ${fq_file1%%.fastq*}_val_1.fq*_bismark_pe.sam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

#bedgraph creation
bismark2bedGraph --CX CpG* -o ${fileID}_CpG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
bismark2bedGraph --CX CHG* -o ${fileID}_CHG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
bismark2bedGraph --CX CHH* -o ${fileID}_CHH.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
cd ../
mkdir 5_output_files
mv 4_bismark_alignment/*.bed* 5_output_files

#100bp window creation

perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CpG* 100 0 ${fileID}_CpG 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CpG*.wig 5_output_files/${fileID}_CpG_100bp.wig
perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CHG* 100 0 ${fileID}_CHG 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CHG*.wig 5_output_files/${fileID}_CHG_100bp.wig
perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CHH* 100 0 ${fileID}_CHH 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CHH*.wig 5_output_files/${fileID}_CHH_100bp.wig

echo "#####################"
echo "providing pipeline metrics to wgbs pipeline logfile..."
echo "#####################"

#get all relevant numbers for final log summary
bismark_version=$(bismark --version | grep "Bismark Version:" | cut -d":" -f2 | tr -d ' ')
samtools_version=$(samtools 3>&1 1>&2 2>&3 | grep "Version:" | cut -d' ' -f2 | tr -d ' ')

map_ef=$(grep 'Mapping efficiency:' 4_bismark_alignment/${fq_file1%%.fastq*}_val_1.fq*_bismark_PE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
unique_aln=$(grep 'Number of paired-end alignments with a unique best hit:' 4_bismark_alignment/${fq_file1%%.fastq*}_val_1.fq*_bismark_PE_report.txt  | cut -d: -f2 | tr -d '\t')
no_aln=$(grep 'Sequence pairs with no alignments under any condition:' 4_bismark_alignment/${fq_file1%%.fastq*}_val_1.fq*_bismark_PE_report.txt  | cut -d: -f2 | tr -d '\t')
multi_aln=$(grep 'Sequence pairs did not map uniquely:' 4_bismark_alignment/${fq_file1%%.fastq*}_val_1.fq*_bismark_PE_report.txt  | cut -d: -f2 | tr -d '\t')
cpg_per=$(grep 'C methylated in CpG context:' 4_bismark_alignment/${fq_file1%%.fastq*}_val_1.fq*_bismark_PE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chg_per=$(grep 'C methylated in CHG context:' 4_bismark_alignment/${fq_file1%%.fastq*}_val_1.fq*_bismark_PE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chh_per=$(grep 'C methylated in CHH context:' 4_bismark_alignment/${fq_file1%%.fastq*}_val_1.fq*_bismark_PE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)



if [[ $fq_file1 == *gz* ]];then
	raw_reads=$(zcat 0_rawfastq/*.gz | wc -l)
	raw_reads=$(($raw_reads / 4 ))
fi

if [[ $fq_file1 != *gz* ]];then
	raw_reads=$(wc -l < 0_rawfastq/$fq_file1)
	raw_reads=$(($raw_reads / 4 ))
fi


if [[ $fq_file1 == *gz* ]];then
	flt_reads=$(zcat 2_trimgalore/*.gz | wc -l)
	flt_reads=$(($flt_reads / 4))
fi

if [[ $fq_file1 != *gz* ]];then
	flt_reads=$(wc -l < 2_trimgalore/${fq_file1%%.fastq*}_val_1.fq*)
	flt_reads=$(($flt_reads / 4))
fi

#add it to the full pipeline logfile
printf "${dow}\t${fq_file1}\t${fileID}\t${genome_path##../}\t${type:1}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n"
printf "${dow}\t${fq_file1}\t${fileID}\t${genome_path##../}\t${type:1}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n" >> $HOME/wgbs_se_pipeline_analysis_record.log

echo "####################"
echo "compressing all sam files..."
echo "####################"
#compress sam and unsorted bam files
find -name "*.sam" | xargs pigz
find -name "*pe.bam" | xargs rm

fi

######################################
# PAIRED END PBAT epignome
######################################

if [ "$1" == "-pese" ];then

if [ "$#" -ne 5 ]; then
echo "Missing required arguments for paired-end pbat!"
echo "USAGE: wgbs_se_pipelinev0.2.sh <-pese> <in fastq R1> <in fastq R2> <path to bismark genome folder> <fileID for output files>"
exit 1
fi
#gather input variables
type=$1; #identifying paired end or single end mode
fq_file1=$2; #the input R1 fastq reads
fq_file2=$3; #R2 reads
genome_path=$4; #the path to the genome to be used (bismark genome prepped)
fileID=$5;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing Bismark PBAT paired-end & single-end alignment with the following parameters:"
echo "Type: $type"
echo "Input File: $fq_file1"
echo "Path to bismark genome folder: $genome_path"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo ""
echo "Full logfile of steps: ${fileID}_logs_${dow}.log"
echo "##################"

#develop directory tree
mkdir ${fileID}_wgbspipeline_${dow}
mv $fq_file1 ${fileID}_wgbspipeline_${dow}
mv $fq_file2 ${fileID}_wgbspipeline_${dow}
cd ${fileID}_wgbspipeline_${dow}

#fastqc
mkdir 1_fastqc
fastqc $fq_file1 $fq_file2 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq_file1%%.fastq*}_fastqc* 1_fastqc #
mv ${fq_file2%%.fastq*}_fastqc* 1_fastqc #


#trim_galore
mkdir 2_trimgalore
cd 2_trimgalore
echo "Trimming first 6bp from reads per suggested protocol for Epignome / BPBAT libraries..."
trim_galore --paired --clip_r1 6 --clip_r2 6 ../$fq_file1 ../$fq_file2 2>&1 | tee -a ../${fileID}_logs_${dow}.log
cd ../

#fastqc_again
mkdir 3_trimmed_fastqc
fastqc 2_trimgalore/${fq_file1%%.fastq*}_val_1.fq* 2_trimgalore/${fq_file2%%.fastq*}_val_2.fq* 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 2_trimgalore/${fq_file1%%.fastq*}_val_1_fastqc* 3_trimmed_fastqc
mv 2_trimgalore/${fq_file2%%.fastq*}_val_2_fastqc* 3_trimmed_fastqc


mkdir 0_rawfastq
mv $fq_file1 0_rawfastq
mv $fq_file2 0_rawfastq

#bismark
mkdir 4_bismark_alignment
cd 4_bismark_alignment
bismark -n 2 -l 20 --bowtie1 --sam --un ../../$genome_path -1 ../2_trimgalore/${fq_file1%%.fastq*}_val_1.fq* -2 ../2_trimgalore/${fq_file2%%.fastq*}_val_2.fq* 2>&1 | tee -a ../${fileID}_logs_${dow}.log

#part 2 - SE directional on unmapped R1
bismark -n 2 -l 20 --bowtie1 --sam ../../$genome_path  ${fq_file1%%.fastq*}_val_1.fq_unmapped_reads_1.fq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

#part 3 - SE directional pbat on unmapped R2
bismark -n 2 -l 20 --bowtie1 --sam  ../../$genome_path  ${fq_file2%%.fastq*}_val_2.fq_unmapped_reads_2.fq 2>&1 | tee -a ../${fileID}_logs_${dow}.log
#bismark -n 2 -l 20 --bowtie1 --sam --pbat ../../$genome_path  ${fq_file2%%.fastq*}_val_2.fq_unmapped_reads_2.fq 2>&1 | tee -a ../${fileID}_logs_${dow}.log

#methylation extraction PE
bismark_methylation_extractor --comprehensive --report --buffer_size 8G -p --no_overlap --gzip ${fq_file1%%.fastq*}_val_1.fq*_bismark_pe.sam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

#methylation extraction SE1
bismark_methylation_extractor --comprehensive --report --buffer_size 8G -s --gzip ${fq_file1%%.fastq*}_val_1.fq*_unmapped_reads_1.fq_bismark.sam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

#methylation extraction SE2
bismark_methylation_extractor --comprehensive --report --buffer_size 8G -s --gzip ${fq_file2%%.fastq*}_val_2.fq*_unmapped_reads_2.fq_bismark.sam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

#make bam files for each line of output:
samtools view -b -S -h ${fq_file1%%.fastq*}_val_1.fq*_bismark_pe.sam > ${fq_file1%%.fastq*}_val_1.fq_bismark_pe.bam
samtools view -b -S -h ${fq_file1%%.fastq*}_val_1.fq*_unmapped_reads_1.fq_bismark.sam > ${fq_file1%%.fastq*}_val_1.fq_bismark_se1.bam
samtools view -b -S -h ${fq_file2%%.fastq*}_val_2.fq*_unmapped_reads_2.fq_bismark.sam > ${fq_file2%%.fastq*}_val_2.fq_bismark_se2.bam

#make combined bam file
samtools merge -h ${fq_file1%%.fastq*}_val_1.fq*_bismark_pe.sam ${fq_file1%%.fastq*}_bismark_combined.bam ${fq_file1%%.fastq*}_val_1.fq_bismark_pe.bam ${fq_file1%%.fastq*}_val_1.fq_bismark_se1.bam ${fq_file2%%.fastq*}_val_2.fq_bismark_se2.bam
#sort combined bam file
samtools sort ${fq_file1%%.fastq*}_bismark_combined.bam ${fq_file1%%.fastq*}_bismark_combined.sorted 2>&1 | tee -a ../${fileID}_logs_${dow}.log
#index combined bam file
samtools index ${fq_file1%%.fastq*}_bismark_combined.sorted.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
#make combined sam file
#samtools view -h ${fq_file1%%.fastq*}_bismark_combined.sorted.bam > ${fq_file1%%.fastq*}_bismark_combined.sorted.sam

#merge the met_extract results from each alignment
zcat CpG*.txt.gz > CpG_context_${fileID}_merged.txt
zcat CHG*.txt.gz > CHG_context_${fileID}_merged.txt
zcat CHH*.txt.gz > CHH_context_${fileID}_merged.txt

#100bp window creation
cd ../
mkdir 5_output_files
perl $HOME/scripts/C_context_window_SREedits.pl 4_bismark_alignment/CpG* 100 0 ${fileID}_CpG 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CpG*.wig 5_output_files/${fileID}_CpG_100bp.wig
bismark2bedGraph --CX CHG* -o ${fileID}_CHG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
bismark2bedGraph --CX CHH* -o ${fileID}_CHH.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
cd ../
mkdir 5_output_files
mv 4_bismark_alignment/*.bed* 5_output_files

#100bp window creation

perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CpG* 100 0 ${fileID}_CpG 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CpG*.wig 5_output_files/${fileID}_CpG_100bp.wig
perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CHG* 100 0 ${fileID}_CHG 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CHG*.wig 5_output_files/${fileID}_CHG_100bp.wig
perl $HOME/NGS-scripts/C_context_window_SREedits.pl 4_bismark_alignment/CHH* 100 0 ${fileID}_CHH 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CHH*.wig 5_output_files/${fileID}_CHH_100bp.wig

echo "#####################"
echo "providing pipeline metrics to wgbs pipeline logfile..."
echo "#####################"

#get all relevant numbers for final log summary
bismark_version=$(bismark --version | grep "Bismark Version:" | cut -d":" -f2 | tr -d ' ')
samtools_version=$(samtools 3>&1 1>&2 2>&3 | grep "Version:" | cut -d' ' -f2 | tr -d ' ')

map_ef=$(grep 'Mapping efficiency:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
unique_aln=$(grep 'Number of alignments with a unique best hit from the different alignments:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
no_aln=$(grep 'Sequences with no alignments under any condition:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
multi_aln=$(grep 'Sequences did not map uniquely:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
cpg_per=$(grep 'C methylated in CpG context:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chg_per=$(grep 'C methylated in CHG context:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chh_per=$(grep 'C methylated in CHH context:' 4_bismark_alignment/${fq_file%%.fastq*}_trimmed.fq*_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)

if [[ $fq_file == *gz* ]];then
	raw_reads=$(zcat 0_rawfastq/*.gz | wc -l)
	raw_reads=$(($raw_reads / 4 ))
fi

if [[ $fq_file != *gz* ]];then
	raw_reads=$(wc -l < 0_rawfastq/$fq_file)
	raw_reads=$(($raw_reads / 4 ))
fi


if [[ $fq_file == *gz* ]];then
	flt_reads=$(zcat 2_trimgalore/*.gz | wc -l)
	flt_reads=$(($flt_reads / 4))


if [[ $fq_file1 == *gz* ]];then
	raw_reads=$(zcat 0_rawfastq/*.gz | wc -l)
	raw_reads=$(($raw_reads / 4 ))
fi

if [[ $fq_file1 != *gz* ]];then
	raw_reads=$(wc -l < 0_rawfastq/$fq_file1)
	raw_reads=$(($raw_reads / 4 ))
fi


if [[ $fq_file1 == *gz* ]];then
	flt_reads=$(zcat 2_trimgalore/*.gz | wc -l)
	flt_reads=$(($flt_reads / 4))
fi

if [[ $fq_file1 != *gz* ]];then
	flt_reads=$(wc -l < 2_trimgalore/${fq_file1%%.fastq*}_val_1.fq*)
	flt_reads=$(($flt_reads / 4))
fi

#add it to the full pipeline logfile
echo "####################"
echo "PE mapping results"
echo "####################"
printf "${dow}\t${fq_file1}\t${fileID}\t${genome_path##../}\t${type:1}_PE\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n"
echo "####################"
echo "Unmapped Read 1 SE mapping results"
echo "####################"
printf "${dow}\t${fq_file1}\t${fileID}\t${genome_path##../}\t${type:1}_SE1\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef2}\t${unique_aln2}\t${no_aln2}\t${multi_aln2}\t${cpg_per2}\t${chg_per2}\t${chh_per2}\n"
echo "####################"
echo "Unmapped Read 2 SE PBAT mapping results"
echo "####################"
printf "${dow}\t${fq_file2}\t${fileID}\t${genome_path##../}\t${type:1}_SE2pbat\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef3}\t${unique_aln3}\t${no_aln3}\t${multi_aln3}\t${cpg_per3}\t${chg_per3}\t${chh_per3}\n"



printf "${dow}\t${fq_file1}\t${fileID}\t${genome_path##../}\t${type:1}_PE\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n" >> $HOME/wgbs_se_pipeline_analysis_record.log
printf "${dow}\t${fq_file1}\t${fileID}\t${genome_path##../}\t${type:1}_SE1\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef2}\t${unique_aln2}\t${no_aln2}\t${multi_aln2}\t${cpg_per2}\t${chg_per2}\t${chh_per2}\n" >> $HOME/wgbs_se_pipeline_analysis_record.log
printf "${dow}\t${fq_file1}\t${fileID}\t${genome_path##../}\t${type:1}_SE2pbat\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef3}\t${unique_aln3}\t${no_aln3}\t${multi_aln3}\t${cpg_per3}\t${chg_per3}\t${chh_per3}\n" >> $HOME/wgbs_se_pipeline_analysis_record.log

echo "####################"
echo "compressing all sam files..."
echo "####################"
#compress sam and unsorted bam files
find -name "*.sam" | xargs pigz
find -name "*_merged.txt" | xargs pigz
fi


#!/bin/bash
set -u
# Based on SRE Bisulfite sequence analysis pipeline
###################
# This script it designed to take a single end fastq file and process it through the 
# bismark aligner call methylated cytosines, and develop per-c bed files and 100bp 
# window wig files for CpG, CHG, and CHH methylation levels
#
# Genome indexing
# Bowtie: bismark_genome_preparation /path/to/genome
# Bowtie2: bismark_genome_preparation --bowtie2 /path/to/genome
