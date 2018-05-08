##### DEG_calling_edgeR #####

# this script contains R code used for calling differential gene expression using edgeR package
# takes count files outputted by featureCounts and calculates differential gene expression between samples/groups of interest
# adapted from pedrocrisp/RNAseq-old and dtrain/NGS-scripts/RNA/edgeR_DEG.r
# need to run script manually
# use read_summarisation_featureCounts.sh to generate count files from RNA-seq data


##### Set up

# run R
R

# Install edgeR (version 3.20.9) if needed
source("http://bioconductor.org/biocLite.R")biocLite(c('edgeR'))

# load edgeR library
library(edgeR)


##### read in files

### get list of featureCount count files that are of interest (i.e. that you want to call DEGs between)
# '^' indicated start of string; '(.*)' indicates anything; '$' indicated end of string
files=list.files(pattern = "^Exp532_(.*)3-r(.*).counts$")


### read in files of interest using readDGE function
# this function reads and merges a number of text files of count data to create a DGList object (object type that is used by edgeR)
#	assumes gene identifies are listed in first column and counts in second column of file
#	creates a table containing rows for genes and columns for samples (genes not found in a sample are given a count of 0)
#	files assumed to be tab delimited and contain column headers
#	if input files are dataframes, file names and group identifiers can be given as columns 'files' and 'groups' 
#	options are as follows:
#		files=vector of file names or a dataframe of sample info containing a 'files' column
#		path=path to directory containing the files of interest (default=current directory)
#		columns=vector containing columns of input files that contain gene names (1) and counts (2)
#		group=vector containing experimental group that each file belongs to (i.e. 0 or 1)
#		labels=vector containing names for each file (default is the file name)
#		can also use other arguments from R's read.delim function (in particular, skip argument to skip the first file of featureCounts file which contains experimental info from featureCounts call, but is not needed for edgeR)

# create 'group' vector containing experimental grouping information (1=control, 2=treated)
groups=c(1,1,1,2,2,2)

# create 'labels' vector containing names for each file - labelled as "genotype"_"replicate" (extracted from file names using strsplit function)
labs <- sapply(strsplit(files, "_"), function(l) paste0(l[2],"_", l[3]))

# read in files and create DGEList object using readDGE function
# for this, columns 1 (geneIDs) and 7 (read counts/gene) from featureCounts output are input
# this command creates an output file containing two columns - samples and counts
#	sample column contains experimental info in sub columns: files, lib.size, group, and norm.factors
#	count column contains count info 

DGEs=readDGE(files, columns=c(1,7), group=groups, labels=labs, skip=1)



##### Filter data

### remove rRNA contamination
# reads present from rRNA contribute to the library size (total count/sequence depth for each library/sample)
#	Since library size is used to make adjustments when calling differentially expressed genes (discussed later on), it is important to remove rRNA prior to DGE calling in order to eliminate effects of an inflated library size
#	most rRNA should be removed from library during poly-A selection, though need to also remove reads from any rRNAs that weren't removed during poly-A selection
# code adapted from dtrain/RNA/edgeR_DEG.r

# get list of genes present in DGEList object (which are listed as the rownames of the DGE$counts column)
geneIDs <- as.character(rownames(DGEs$counts))

# read in list of Arabidopsis rRNA IDs (this is generated using GTF_to_annotation.sh script and should be located in annotations directory)
rRNA <- read.delim("~/annotations/Ara11_rRNA.txt", head=F)
rRNA <- as.character(rRNA$V1)

# identify rRNA reads present in DGEs library using match command (matches entries in one list to those in another)
rRNA_contam <- match(rRNA, geneIDs)


# count number of rRNA reads identified in DGEs library by subsetting DGE$counts column to rRNAs
rRNA_counts <- DGEs$counts[rRNA_contam, ]

# count the number of total rRNA reads in the DGE library by summing the counts for each individual rRNA 
rRNA_sum <- colSums(rRNA_counts)


# calculate the rate of rRNA contamination in the DGE library (this is calculated as total rRNA count/DGE library size)
rRNA_rates <- (rRNA_sum/DGEs$samples$lib.size)

# plot rRNA contamination rates
pdf("20180428_Exp532_Col-0_WLRSII_rRNA_contamination", paper="a4r")
plot(rRNA_rates, main = "rRNA contamination in libraries", ylab = "rRNA abundance (% of total mapped reads)",
     xlab = "Sample")
dev.off()

# remove rRNA reads from library
DGEs$counts <- DGEs$counts[-rRNA_contam, ]


### remove chloroplastic and mitochondrial transcript contamination
# these should not be present in mRNA-seq data (as the plastids lack capacity for polyadenylation) and so are considered as contamination

# get list of genesIDs in count table (after rRNA removal)
geneIDs <- as.character(rownames(DGEs$counts))

# get annotation file containing mitochondrial and chloroplastic geneIDs
anno=read.delim("/home/bethanys/annotations/Araport11_genes.sorted.bed", head=F)

# subset annotation to genes from chloroplast and mitochondria
p_genes=subset(anno, anno$V1 == "ChrM" | anno$V1 == "ChrC")

# get list of plastid genes as character vector
p_genes=as.character(p_genes$V4)

# match these to geneID list to get vector numbers of counts to be removed
plastid_contam=match(p_genes, geneIDs)

# remove NA values
plastid_contam=na.omit(plastid_contam)

# remove plastid reads from library
DGEs$counts <- DGEs$counts[-plastid_contam, ]


### remove lowly expressed genes 
# genes with low counts across libraries give little evidence for differential expression (i.e. a gene must be expressed at some min level for it to be likely to be translated/be biologically relevant)
# pronounced discreteness of these genes also interferes with statistics used for DGE analysis
# genes should be filtered out if they can't be expressed in all samples for any of the conditions
# typically, genes required to have count of 5-10 to be considered
# filter based on CPM not counts directory (need to account for differences in library size between samples)

# sum counts over lowly expressed genes
# row sums condenses rows of a DGEList object and sums counts over specific groups of genes
# arguments are as follows
#		x = DGEList object
#		group = vector or factor level giving the grouping
#		reorder = if TRUE, rownames in output DGEList are sorted as sort(unique(group)); if FALSE, they are given in order that groups were encountered
#		na.rm = TRUE or FALSE whether NA values should be discarded?
# here, genes are kept if they have a CPM greater than 1 in at least 3 samples

keep <- rowSums(cpm(DGEs) > 1) >= 3

# subset DGE library to genes to be kept 
DGEs <- DGEs[keep, ]


### recalculate library size post-filtering
# sequencing depth of mRNA samples of interest affects the read counts for genes, and so needs to be adjusted for when calling differential expression as different libraries (samples) likely have different sequencing depths
# edgeR takes this account automatically during modelling of fold-changes or p-value calculations
# library size thus needs to be re-calculated following the above filtering to improve accuracy of DGE analysis (although changes are unlikely to be very large...)

# re-calculate library sizes based on kept transcripts
DGEs$samples$lib.size <- colSums(DGEs$counts)



##### Normalisation
# edgeR only looks at relative changes in expression between groups/conditions rather than quantification of expression levels
#	this means that technical factors unrelated to experimental condition should not influence differential expression analysis
#	thus, edgeR normalisation is only needed for technical factors that have sample-specific effects
#	technical factors that may need to be normalised include sequencing depth, RNA composition, and to a lesser extent CG content and gene length
# sequencing depth is accounted for by edgeR automatically (see above)
# RNA composition - RNA-seq gives a measure of relative abundance of each gene but doesnt provide a measure of total RNA output per cell
#	problem when small numbers of genes are highly expressed in one sample but not in another (as they can take up a high proportion of total library size and so cause the rest of the genes in the sample/library to be under-sampled and so may appear as falsely down-regulated)

### normalise library based on RNA composition
# calcNormFactors function normalises RNA composition using a set of scaling factors for the library sizes that minimise the log-fold changes between samples for the most genes
#	it calculates normalization factors to scale the raw library sizes by
# default method for calculating scaling factors is trimmed-mean M-values (TMM) between each pair of samples
# an effective library size is calculated based on the scaling factors, which is used for DGE calling
# TMM normalisation is recommended for most RNA-seq analyses (when most genes are believed to be not differentially expressed between any pair of the samples

DGEs_TMM <- calcNormFactors(DGEs, method = "TMM")


### estimating dispersion 
# NEED TO FILL IN THIS INFO

# generate design matrix
design <- model.matrix(~groups)

# estimate common, trended, and tagwise dispersion
DGEs_TMM_disp <- estimateDisp(DGEs_TMM, design=design)


# Plot biological coefficient of variation
plotBCV(DGEs_TMM_disp)
groups <- unique(as.character(samplegroups))
n.reps = length(sampleGroups)/length(groups)

# plot multidimensional scaling plot
plotMDS(dge.tmm.disp, dim.plot = c(1,2), col = rep(rainbow(length(groups)), each = n.reps))
plotMDS(dge.tmm.disp, dim.plot = c(2,3), col = rep(rainbow(length(groups)), each = n.reps)) 



##### test for differential gene expression

### single factor exact tests
# once negative binomial model is fitted/dispersion estimated, differential expression can be analysed using an exact test
# exact test computes exact p-values by summing all sums of counts that have a "probability less than the probability under the null hypothesis of the observed sum of counts"
# exact tests only applicable to experiments with a single factor
# can use common dispersion or tagwise dispersion approaches 
# exactTest() function computes genewise exact tests for differences in the means between two groups of negative-binomially distributed counts
# arguments are as follows
#	object = DGEList class object
#	disperson = numeric vector of dispersion values or character string indicating to take dispersions from the object ("common", "trended", "tagwise"). Default is "auto" which uses the most complex dispersion found in the object
#	pairs = numeric or character vector providing the pair of groups to be compared (if default, the first two levels of object$samples$group will be used.  The first group listed in the pair is used as the baseline for the comparison (i.e. A-B => +ve fold change is upregulated in group B compared to A)
# trended dispersions are generally used to account for empirical mean-variance relationships

ET <- exactTest(DGEs_TMM_disp, dispersion="trended")


### find top DEGs from exact test output
# topTags function extracts the top DEGs from the exact test output dataframe for a given pair of groups
# DEGs can be ranked by  p-value or absolute log-fold change
# arguments are as follows:
#	object = DGEExtract object outputted from exactTest function
#	n = number of tags to display
#	adjust.method = method used for adjusting p-values for multiple testing
#	sort.by = whether to sort top tags by p-value or absolute log-fld change, or none ("PValue", "logFC", or "none")
#	p.value = max value for adjusted p-values (i.e. tags with higher p-values are excluded)
# p-value adjustment methods could be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" or "none"
#	here "fdr" is used which controls the false discovery rate (expected proportion of false discoveries amongst the rejected hypothesis)

ET_TTs <- topTags(ET, adjust.method = "fdr", sort.by="logFC", p.value=0.05)


### Find genes that are significantly differentially expressed
# decideTestsDGE() function uses multiple testing procedure and significance level cutoff to find significantly differentially expressed genes in a genewise tests object (i.e. a DGEExact object)
# for expression to be considered as significantly different, the gene must have an adjusted p-value below the specified p-value and a log2-FC greater than the specified log2-fold change value
# arguments are as follows:
#	object = DGEExact object containing p-values and log-FCs
#	adjust.method = p-value adjustment method ("none", "BH"/"fdr", "BY" and "holm")
# 	p.value = number between 0 and 1 giving the required error rate/false discovery rate.
#	lfc = number specifying the min absolute log2-fold-change required
# -1 = significantly down-regulated, 0 = not significant, 1 = significantly up-regulated

DGEs_sig <- decideTestsDGE(ET, adjust.method = "fdr", p.value = 0.05)

# generate summary of significant DEGs (tells you how many genes are up or down-regulated)
summary(DGEs_sig)



##### get rpkm values 

### use Araport11 annotation file to obtain gene lengths

library(tidyverse)
ara11 <- read.delim("/home/bethanys/annotations/Araport11_genes.sorted.bed", head=F) %>%
	mutate(length = V3 - V2) %>%
	filter(V4 %in% geneIDs)
gene_lengths <- ara11$length[ara11$V4 %in% rownames(DGEs_TMM_disp$counts)]

### calculate rpkm for each gene per sample
rpkms <- rpkm(DGEs_TMM_disp, gene.length=gene_lengths, dispersion=DGEs_TMM_disp$trended.dispersion)

### calculate average rpkm for each gene per group
rpkms_groups <- rpkmByGroup(DGEs_TMM_disp, gene.length=gene_lengths, dispersion=DGEs_TMM_disp$trended.dispersion)



##### generate table of DGEs for outputting

### column bind significant DGEs to exact test table output
ET_full=cbind(ET$table, DGEs_sig)


### add rpkm values to ET_full dataframe
ET_full=cbind(ET_full, rpkms)
ET_full=cbind(ET_full, rpkms_groups)


### subset ET_full to DGEs and write out file
ET.out=subset(ET_full, ET_full$'1+2' != 0)

### write out file
write.table(ET.out, "Exp532_WLRSII_Col-0_vs_strs2_trt_DGEs", col.names=TRUE, row.names=TRUE, quote=FALSE, sep='\t')

