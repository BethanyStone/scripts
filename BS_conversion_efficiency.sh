#!/bin/bash
set -e
set -u

##### BISULFITE CONVERSION EFFICIENCY

# based on dtrain/conversion_rate_check.sh
# this script checks bisulfite conversion efficiency by calculating the percentage of unconverted cytosines in the chloroplastic and mitochondrial genome - these genomes should be unmethylated, so perfect conversion should be 0% and greater than this indicates over conversion.
# uses bismark coverage file as input (chromosome start_position end_position %methylation count_methylated count_unmethylated
# % methylated calculated as:
#	(# methylated cytosines)/(# total cytosines) * 100
# % conversion efficiency:
#	 100 - ((# methylated cytosines)/(# total cytosines) * 100)
# this is done for chloroplastic and mitochondrial genomes

### usage
if [ "$#" -lt 1 ]; then
echo "Missing required arguments!"
echo "USAGE: BS_conversion_efficiency.sh bismark_CHH_cov_file"
exit 1
fi


### define input variables
cov_file=$1

### unzip cov_file
if [[ ${cov_file} == *gz ]];then
	gunzup ${cov_file}
fi


### check chromosome naming in input coverage file
# only need to run this if chromosome naming is unknown... uncomment if needed
#awk -F '\t' '{print $1}' ${cov_file} | sort | uniq -c | sort -nr


### calculate conversion efficiecies
echo "%met_chrC"
grep "chloroplast"  ${cov_file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print (met / total)*100}'
echo "chrC_efficiency"
grep "chloroplast"  ${cov_file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
echo "%met_chrM"
grep "mitochondria"  ${cov_file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print (met / total)*100}'
echo "chrM_efficiency"
grep "mitochondria"  ${cov_file} | awk '{ met+= $5} { unmet += $6} { total = met + unmet } END {print 100-((met / total)*100)}'
echo "done"


### gzip cov_file
if [[ ${cov_file} != *gz ]];then
	gzip ${cov_file}
fi

