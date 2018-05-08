##### GO_analysis_topGO.R #####

# quick script to perform simple gene ontology analysis for gene sets of interest (i.e. DEGs or DAGs)
# don't fully understand workings of topGO, but following code based on this handy post: https://www.biostars.org/p/250927/
# https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
# https://bioconductor.org/packages/release/bioc/manuals/topGO/man/topGO.pdf
# script needs to be altered/run manually...


##### install packages, run R and load topGO library

### package installation
# to install topGO
source("https://bioconductor.org/biocLite.R")
biocLite("topGO")

# to install biomart for R
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

# to install Rgraphviz
source("https://bioconductor.org/biocLite.R")
biocLite("Rgraphviz")


### run R
R


### load packages
# load toGO library
library(topGO)

# load biomart library
library(biomaRt)



##### create a geneID2GO annotation 

# the geneID2GO object contains a list of which gene IDs map to which GO IDs
# genes IDs need to match the formatting of the gene IDs used in the DEG/DMR analysis (i.e. those in the list of genes you're interested in...)
# should be given as a named list of characters vectors
#	list names give the gene IDs
#	list elements give a character vector containing the GO identifiers for each listed gene 
# ensembl gene annotations are retrieved using biomaRt package


### get ensembl annotation information using biomaRt
# to check which annotations are available
biomaRt::listMarts(host="plants.ensembl.org")
datasets=biomaRt::listDatasets(biomaRt::useMart(biomart="plants_mart",host="plants.ensembl.org"))
table(datasets$dataset)

# load in the ensembl athaliana gene annotation from biomart using useMart function
anno <- biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')

# retrieve geneIDs and GO IDs using getBM function and lists in a tab-delimited table
# this is the main query function in biomaRt, with arguments as follows:
#	attribute = vector of attributes to retrieve
#	mart = object of mart class created using useMart (in this case, "anno")
G2GO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = anno)


### format resulting G2GO dataframe to create geneID2GO object
# remove any geneIDs that have missing GO terms
G2GO <- G2GO[G2GO$go_id != '',]

# create the final geneID2GO annotation object by changing table to list 
# note - entries need to be as characters
geneID2GO <- by(G2GO$go_id, G2GO$ensembl_gene_id, function(x) as.character(x))

# check structure of geneID2GO object is correct
head(geneID2GO)



##### create gene lists for GO analysis

### create a list of all gene IDs (using geneIDs from G2GO object)
geneNames=sort(unique(as.character(G2GO$ensembl_gene_id)))

### create a list of genes of interest (i.e. DEGs or DAGs)
# read in your file containing genes of interest
DEGs=read.delim("Exp532_WLRSII_Col-0_ctrl_vs_trt_DGEs", head=T, row.names=NULL)

# relabel column containing geneIDs (if needed)
colnames(DEGs)[1]="ID"

# create character vector containing geneIDs of genes of interest (GOIs)
GOIs=as.character(DEGs$ID)

### generate geneList object
# this contains integers (0 or 1) for all genes indicating whether genes 'interesting' or not
# 0 = not GOI, 1 = GOI
geneList <- factor(as.integer(geneNames %in% GOIs))

# add geneIDs as names to vector of GOI intergers
names(geneList) = geneNames


##### create a topGOdata object
# contains all necessary info for GO analysis
#	gene list, list of interesting genes, gene scores (if available), GO ontology information
# to build an object you need:
#	list of geneIDs (and gene-wise scores, i.e. p value for diff expresssion/diff methylation)
# 	mapping between gene geneIDs and GO terms (need to specify name of the annotation being used)
#	GO hierachical structure obtained from the GO.db package
# don't fully get what this is doing... but it works!

GOdata <- new("topGOdata", ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)



##### Perform tests for GO enrichment in GOIs

# here an enrichment test based on gene counts is used - this only requires a list of interesting genes
### use runTest function to run GO enrichment test
# runTest is a wrapper for getSigGroups function
# arguments are as follows
#	object = topGOdata object
#	test.stat = defines the test statistic
#	algorithm = specifies what algorithm to use for test ("classic", "elim", "weight", "weight01", "lea", "parentchild")
# here, a classic algorithm and fisher is used to test for over-representation of GO terms within the group of GOIs (as per https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf)

rslts <- runTest(GOdata, algorithm = "classic", statistic = "fisher")


### generate a table of GO enrichment results

# use GenTable function to generate a table of GO enrichment results (by GO term, not by gene... there's another function to do this if needed)
# arguments are as follows
#	object = topGOdata class object
#	ranksOf = scores ordered by ranks of the specified result
#	topNodes = number of GO terms to be included in the table
#	numChar = GO term descriptions truncated to show this number of characters

rsltsTab <- GenTable(object = GOdata, ranksOf = "classic", Classic = rslts, topNodes = 20, numChar=200)

# write this out as a table
write.table(rsltsTab, "Exp532_WLRSII_Col-0_ctrl_vs_trt_DGEs_GO.txt", sep='\t', row.names=F, col.names=T, quote=F)


### generate a GO graph of enriched GO terms

# use printGraph function, with arguments as follows:
#	object = topGOdata class object
#	result = topGOresult class object
#	firstSigNodes = number of top scoring nodes
#	fn.prefix = file name prefix
#	useInfo = all (use all nodes)
#	pdfSW = TRUE (prints as PDF file)
# significant nodes are represented as rectangles in the graph	

printGraph(GOdata, rslts, firstSigNodes = 10, fn.prefix = "Exp532_WLRSII_Col-0_ctrl_vs_trt_DGEs", useInfo = "def", pdfSW=TRUE)
