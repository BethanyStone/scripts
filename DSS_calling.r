# Script to perform DSS DMR calling - based on dtrain DSS_calling.r
# Need to run R code manually
# Files of interest need to be transformed into appropriate format for DSS using DSS_file_prep.r

# install DSS (if needed)
# source("http://bioconductor.org/biocLite.R")
# biocLite("DSS")

# start R
R

# Define following
context = "ENTER CONTEXT YOU'RE INTERESTED IN HERE"
pvalue = 0.05
delta = "ENTER DELTA VALUE" # 0.4 for CpG, 0.2 for CHG, 0.1 for CHH

library(DSS)

# Read in correctly formatted files
files <- list.files(pattern = paste0(context,".output"))

# Define sample groups
group1 <- files[c(1,2,3)]
group2 <- files[c(4,5,6)]

# read input files in DSS format (chr, pos, N, X)
dat1.1 <- read.delim(unlist(group1)[1])
dat1.2 <- read.delim(unlist(group1)[2])
dat1.3 <- read.delim(unlist(group1)[3])

dat2.1 <- read.delim(unlist(group2)[1])
dat2.2 <- read.delim(unlist(group2)[2])
dat2.3 <- read.delim(unlist(group2)[3])

# setup bsseq object
BSobj <- makeBSseqData(list(dat1.1,dat1.2,dat1.3,dat2.1,dat2.2,dat2.3),sampleNames=c("C1","C2","C3","N1","N2","N3"))

# Estimation of methylation means with smoothing by moving averages and smaller smoothing window
# uses 200bp spanning windows for smoothing
dmlTest <- DMLtest(BSobj,group1=c("C1","C2","C3"), group2=c("N1","N2","N3"),smoothing=TRUE,smoothing.span=200)

# identify DMLs and write out
dmls <- callDML(dmlTest, delta=delta, p.threshold=pvalue)

# Write dmls to table - need to enter file name of choice
write.table(dmls, "FILE_NAME", quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# identify DMRs based on dmltesting and write out to file
dmrs <- callDMR(dmlTest, delta=delta, minlen=50, minCG=3, pct.sig=0.5, dis.merge=50, p.threshold=pvalue)

#write dmrs to file - need to enter file name manually
write.table(dmrs,"FILE_NAME", quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
