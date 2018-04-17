# DMR ANALYSIS
# bits of code used for dmr analysis following DMR calling with DSS
# not formatted yet, just keeping track of commands used...

##### Hypo vs Hyper DMRs #####
### to count number of hyper (-ve diff) and hypo (+ve diff) DMRs:

# start R
R

# read in DSS output txt files
CpG_dmrs=read.delim("Col-0_vs_strs2_CpG_output.txt_DMRs_CpG_p=0.05.bed", head=T)
CHG_dmrs=read.delim("Col-0_vs_strs2_CHG_output.txt_DMRs_CHG_p=0.05.bed", head=T)
CHH_dmrs=read.delim("Col-0_vs_strs2_CHH_output.txt_DMRs_CHH_p=0.05.bed", head=T)

#count number of rows (equal to number of DMRs)
nrow(CpG_dmrs)
nrow(CHG_dmrs)
nrow(CHH_dmrs)

# count number of DRMS based on sign of diff.methyl column 
# Col-0 = control, strs2 = "trt" - therefore, +ve diff.methyl value = hypoDMR (col>strs2) and -ve diff.methyl value = hyperDMR (strs2>col)
table(sign(CpG_dmrs$diff.Methy))
table(sign(CHG_dmrs$diff.Methy))
table(sign(CHH_dmrs$diff.Methy))


##### DMRS at genes vs TEs #####
### find DMRs intesecting genes or TEs - NOT IN R
# skip first header line
tail -n +2 Col-0_vs_strs2_CpG_output.txt_DMRs_CpG_p=0.05.bed > Col-0_vs_strs2_CpG_output.txt_DMRs_CpG_p=0.05.noheader.bed
tail -n +2 Col-0_vs_strs2_CHG_output.txt_DMRs_CHG_p=0.05.bed > Col-0_vs_strs2_CHG_output.txt_DMRs_CHG_p=0.05.noheader.bed
tail -n +2 Col-0_vs_strs2_CHH_output.txt_DMRs_CHH_p=0.05.bed > Col-0_vs_strs2_CHH_output.txt_DMRs_CHH_p=0.05.noheader.bed

# sort input bed file using sortBed prior to intersect command
sortBed -i Col-0_vs_strs2_CpG_output.txt_DMRs_CpG_p=0.05.noheader.bed > Col-0_vs_strs2_CpG_output.txt_DMRs_CpG_p=0.05.noheader.bed.sorted
sortBed -i Col-0_vs_strs2_CHG_output.txt_DMRs_CHG_p=0.05.noheader.bed > Col-0_vs_strs2_CHG_output.txt_DMRs_CHG_p=0.05.noheader.bed.sorted
sortBed -i Col-0_vs_strs2_CHH_output.txt_DMRs_CHH_p=0.05.noheader.bed > Col-0_vs_strs2_CHH_output.txt_DMRs_CHH_p=0.05.noheader.bed.sorted

# intersect sorted bed files with sorted bed gene and TE annotation files
bedtools intersect -a Col-0_vs_strs2_CpG_output.txt_DMRs_CpG_p\=0.05.noheader.bed.sorted -b ~/annotations/Araport11_genes.sorted.bed ~/annotations/Araport11_TE.sorted.bed -wa -wb -sorted > CpG_genes_vs_TEs_intersect.bed
bedtools intersect -a Col-0_vs_strs2_CHG_output.txt_DMRs_CHG_p\=0.05.noheader.bed.sorted -b ~/annotations/Araport11_genes.sorted.bed ~/annotations/Araport11_TE.sorted.bed -wa -wb -sorted > CHG_genes_vs_TEs_intersect.bed
bedtools intersect -a Col-0_vs_strs2_CHH_output.txt_DMRs_CHH_p\=0.05.noheader.bed.sorted -b ~/annotations/Araport11_genes.sorted.bed ~/annotations/Araport11_TE.sorted.bed -wa -wb -sorted > CHH_genes_vs_TEs_intersect.bed


### to count number of DMRs intersecting genes (1) or TEs (2)
# start R
R

# read in intersect files
CpG=read.delim("CpG_genes_vs_TEs_intersect.bed", head=F)
CHG=read.delim("CHG_genes_vs_TEs_intersect.bed", head=F)
CHH=read.delim("CHH_genes_vs_TEs_intersect.bed", head=F)

# subset to entries intersecting genes (1) or TEs (2)
CpG_genes=subset(CpG, CpG$V10 == 1)
CpG_TEs=subset(CpG, CpG$V10 == 2)

# count hyper/hypo DMRs per feature type
table(sign(CpG_genes$V8))
table(sign(CpG_TEs$V8))

CHG_genes=subset(CHG, CHG$V10 == 1)
CHG_TEs=subset(CHG, CHG$V10 == 2)

table(sign(CHG_genes$V8))
table(sign(CHG_TEs$V8))

CHH_genes=subset(CHH, CHH$V10 == 1)
CHH_TEs=subset(CHH, CHH$V10 == 2)

table(sign(CHH_genes$V8))
table(sign(CHH_TEs$V8))

