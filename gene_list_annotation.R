##### gene_list_annotation.R #####

# quick script to generate arabidopsis gene annotations (including gene descriptions) from ENSEMBL TAIR10.39 annotation and annotate gene lists
# code for parsing gff files from https://support.bioconductor.org/p/24657/
# need to run manually - can alter script depending on what you want in the annotation



##### generate TAIR10 gene annotation

### get ENSEMBL TAIR10 gff annotation file (if you haven't already got it...)

wget ftp://ftp.ensemblgenomes.org/pub/release-39/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.39.gff3.gz


### run R

R


### load in functions for parsing gff/gtf files

# load get attribute fields function
getAttributeField <- function (x, field, attrsep = ";") {
	     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
		             a = strsplit(atts, split = "=", fixed = TRUE)
			              m = match(field, sapply(a, "[", 1))
			              if (!is.na(m)) {
					                   rv = a[[m]][2]
				               }
				               else {
						                    rv = as.character(NA)
					                }
				               return(rv)
					            })
}

# load read gff function
gffRead <- function(gffFile, nrows = -1) {
	     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
		           header=FALSE, comment.char="#", nrows = nrows,
			        colClasses=c("character", "character", "character", "integer",
					     "integer",
					          "character", "character", "character", "character"))
          colnames(gff) = c("seqname", "source", "feature", "start", "end",
			                 "score", "strand", "frame", "attributes")
             cat("found", nrow(gff), "rows with classes:",
		         paste(sapply(gff, class), collapse=", "), "\n")
	       stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
	          return(gff)
}


### read in TAIR10 annotation file

anno=gffRead('Arabidopsis_thaliana.TAIR10.39.gff3')


# check naming of chromosomes (need these to match other annotation files/reference genome file)
table(anno$seqname)

# change naming of chromosomes to be consistent with other files if needed
anno$seqname=gsub("1", "Chr1", anno$seqname)
anno$seqname=gsub("2", "Chr2", anno$seqname)
anno$seqname=gsub("3", "Chr3", anno$seqname)
anno$seqname=gsub("4", "Chr4", anno$seqname)
anno$seqname=gsub("5", "Chr5", anno$seqname)
anno$seqname=gsub("Mt", "mitochondria", anno$seqname)
anno$seqname=gsub("Pt", "chloroplast", anno$seqname)


### subset to annotation feature of interest (in this case, genes)

genes=subset(anno, anno$feature=="gene")


### get attribute information of interest (from attributes column) and add as new column to dataframe

# get geneID
genes$ID=getAttributeField(genes$attributes, 'gene_id')

# get gene name
genes$Name=getAttributeField(genes$attributes, 'Name')

# get gene biotype (i.e. protein coding etc)
genes$Biotype=getAttributeField(genes$attributes, 'biotype')

# get gene description
genes$Description=getAttributeField(genes$attributes, 'description')


### format annotation dataframe

# subset to columns of interest
genes.out=genes[,c('seqname','start','end','ID','score','strand','Name','Biotype','Description')]

# change column names
colnames(genes.out)=c("Chr","Start","End","ID","Score","Strand","Name","Biotype","Description")


### write out annotation file
write.table(genes.out,'TAIR10.39_genes.bed',sep='\t',row.names=F,col.names=F,quote=F)



##### annotating gene lists

### read in list of genes to be annotated (i.e. DEGs, DAGs)

# read in file
DEGs=read.delim("Exp532_WLRSII_Col-0_trt_vs_ctrl_DGEs", head=T, row.names=NULL)

# change any column names if needed
colnames(DEGs)[1]="ID"


### read in annotation file (generated as above)

anno=read.delim("/home/bethanys/annotations/TAIR10.39_genes.bed", head=F)


### merge in annotation info into DEG list

# change columns for merging as necessary
DEGs_annotated=merge(x=DEGs, y=anno, by.x="ID", by.y="V4", all.x=TRUE)


### write out file
write.table(DEGs_annotated, "Exp532_WLRSII_Col-0_trt_vs_ctrl_DGEs_annotated", sep='\t', row.names=F, col.names=T, quote=F)
