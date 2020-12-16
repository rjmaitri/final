#####################################
#'
#' RNA-SEQ quant by Rsubread
#' 12/8/20
#' Final Project Data Wrangling

#
library(Rsubread)
library(edgeR)
library(tidyverse)


####The purpose of the workflow is to find differentially expressed genes from a RNA-seq data of 
# Estrogen-receptor+ breast cancer cell lines (MCF-7), one group treated with 
#fulvestrant (estrogen inhibitor) and the other

#get a list of the raw reads in the directory /rnaseq
fastq.files <- list.files(path = "rnaseq/", pattern = ".fastq$", full.names = TRUE)
fastq.files
##abbreviate sample names
fastq.abbrv <- substr(fastq.files, start = 8, stop = 9)


###VEH_fastq.files was created to resume index alignment after a blue screen of death crash
#Veh_fastq.files <- list.files(path = "rnaseq/", pattern = "^V", full.names = TRUE)
#Veh_fastq.files

#build the index -
#reference genome from ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
# buildindex(basename="Sapiens_Index",reference="GRCh38_latest_genomic.fna", gappedIndex = TRUE, indexSplit = TRUE)

#build index with EBI 
#oaded from the GENCODE database via the following links:
#ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.primary_
#assembly.genome.fa.gz
# buildindex(basename="Ebi_Homo_Sapien_Index",reference="GRCh38.primary_assembly.genome.fa.gz", gappedIndex = TRUE, indexSplit = TRUE)

###############The initial index using ncbi produced very low gene counts - 
#I found out the Ebi index was recommended by the Rsubread authors - 
#Hence the second index build
#the first alignment with EBi index crashed on the 
#V1, the second align() boded below here is to pick up at V1

#maps reads in dataset to reference index
#align.stat <- align(index = "Ebi_Homo_Sapien_Index", readfile1 = fastq.files)

#picking up at V1
#align.stat <- align(index = "Ebi_Homo_Sapien_Index", readfile1 = Veh_fastq.files)



#get annotation file to count mapped reads against
# ann <- getInBuiltAnnotation(annotation = "hg38")

#bild a vector of mapped reads
bam.files <- list.files(path = "rnaseq/", pattern = ".BAM$", full.names = TRUE)


#featureCounts generates a matrix of read counts to genes in each sample
#strandSpecific = 2 because 'TruSeq Stranded mRNAseq Sample Prep kit' (Illumina). 
#Read 1 aligns to the ANTISENSE strand and Read 2 aligns to the SENSE strand,
#this maps the mrna sequence to the annotated exons
fcRS <- featureCounts(bam.files, annot.inbuilt="hg38", strandSpecific = 2)
#tidy sample names(needed for MDS plot labels)
fcRS$targets <- substr(fcRS$targets, start = 1, stop = 2)
####write counts to table (needed for mean/variance dependence)
fcRS_counts <- fcRS$counts
#stranded (SUSPICIOUSLY LOW EXPRESSION)
#fcS <- featureCounts(bam.files, annot.inbuilt="hg38", strandSpecific = 1)
#dimensions - rows are genes, columns are samples
dim(fcRS$counts)
#six samples with 28395 genes

#This outputs the mapping status of the reads
#fc$stat

#SummarizedExperiment HAVENT USED THIS, MIGHT BE USEFUL IF I GET IT TO WORK
#matrix features are accessed with an R function of the same name: assay (counts), rowRanges(geneID) and colData(samples/targets).
##se <- SummarizedExperiment(fc$counts, rowRanges = fc$annotation$GeneID, colData = fc$targets)

####write to file
####JUST USE THIS
#fc_counts <- fcS$counts

#output a txt file of count data
write.table(x=data.frame(fcRS$annotation[,c("GeneID","Length")],fcRS$counts,stringsAsFactors=FALSE),file="RScountsEBI.txt",quote=FALSE,sep="\t",row.names=FALSE)

#For edgeR quasi-likelihood DGE
write.table(x=data.frame(fcRS$annotation[,c("GeneID")],fcRS$counts,stringsAsFactors=FALSE),file="counts.txt",quote=FALSE,sep="\t",row.names=FALSE)

#load data into a table (for count distribution)
EBI_RStrandedCountTable <- read.table(file = 'RScountsEBI.txt', header = TRUE)

#load data into a table for quasi-likelihood DGE
CountTable <- read.table(file = 'counts.txt', header = TRUE)


#######HISTOGRAM############### FUCK YEAHHHH complete facet wrap???
ggplot(EBI_RStrandedCountTable) +
  geom_histogram(aes(x = F1_CAAGCTAG.ACATAGCG_L002_R1_001..1..fastq.subread.BAM), stat = "bin", bins = 250, colour = "lightblue", fill = "orange") +
  xlab("Raw expression counts") +
  ylab("Number of genes")+
  xlim(0,20000)+
  ylim(0,1760)+
  #facet_wrap()+
  theme_dark()
#F1 COUNT DISTRIBUTION(gamma) = total gene amount vs. expression counts shows a typical single cell sequencing
#distribution, with alot of lowly expressed genes and few genes with high expression

#extract null counts as percentage
prop.null <- apply(EBI_RStrandedCountTable, 2, function(x) 100*mean(x==0))
View(prop.null)
print(head(prop.null))
#plot null counts
barplot(prop.null, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, ylab='Samples', xlab='% of null counts')

#Determine which genes have sufficiently large counts to be retained in a statistical analysis.
#default values of filtByExpr have been used for low expression - 
#count.table <- filterByExpr(fc$counts, group = fastq.files)
#fc <- fc$counts[count.table,]

#alot of misses
#compute vectors for mean and variance of fulvestrant-treated samples
mean <- apply(fcRS_counts[,1:3], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance <- apply(fcRS_counts[,1:3], 1, var)
MeanVariancedf <- data.frame(mean, variance)

#scatter-plot to display mean-variance relationship
ggplot(MeanVariancedf) +
  geom_point(aes(x=mean, y=variance), colour = "orange", alpha = .8) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="green")+
  theme_dark()

####looks typical of single cell seq data, many genes with low-level of expression and a few
#high-level expression

#############RSUBUSERS RNA-SEQ DGE ANALYSIS

DGEfulvestrant <- DGEList(counts=fcRS$counts, genes=fcRS$annotation[,c("GeneID","Length")])

####Filter to keep genes which had  more than 10 reads per million mapped
#reads in at least two libraries. **credit Rsubread Users Guide
FilteredGenes <- rowSums(cpm(DGEfulvestrant) > 10) >= 2
#write to dataframe
DGEfulvestrant <- DGEfulvestrant[FilteredGenes,]
##Design matrix. Create a design matrix:

samples <- factor(fcRS$targets)
designmatrix <- model.matrix(~0+samples)
colnames(designmatrix) <- levels(samples)
##Normalization. Perform voom normalization: DOESNT WORK
#what is the math behind normalization in this function????
Normalization <- voom(DGEfulvestrant,designmatrix,plot=T)



#####MULTI-DIMENSIONAL SCALING PLOT
####dimension analysis - change names of samples
#distances on the plot approximate the typical 
#log2 fold changes between the samples.
#This function is a variation on the usual multdimensional 
#scaling (or principle coordinate) plot, in that a distance 
#measure particularly appropriate for the microarray context 
#is used. The distance between each pair of samples 
#(columns) is the root-mean-square deviation (Euclidean distance) 
#for the top top genes. Distances on the plot can be 
#interpreted as leading log2-fold-change, meaning the 
#typical (root-mean-square) log2-fold-change between the samples for the genes that distinguish those samples.
plotMDS(Normalization, labels = samples, xlim=c(-1,1),ylim=c(-1.0,.5))

# Fancy Shmancy multi-dimensional plot i can try later
par(mfrow=c(1,2))
col.cell <- c("purple","orange")[sampleinfo$CellType]
col.status <- c("blue","red","dark green")[sampleinfo$Status]
plotMDS(dgeObj,col=col.cell)
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
title("Cell type")
plotMDS(dgeObj,col=col.status)
legend("topleft",fill=c("blue","red","dark green"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

#The distance between each pair of samples in the MDS plot is 
#calculated as the leading fold change, defined as the root-mean-square 
#of the largest 500 log2-fold changes between that pair of samples.

x <- read.delim("counts.txt",row.names="fcRS.annotation...c..GeneID...")
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

#Fit a glm quasi-likelihood
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)



########results
p-values/heatmap/dge list/gene ontology/pattern pathways

#'Variance-mean dependence was estimated from count tables and 
#'tested for differential expression based on a negative binomial 
#'distribution, using DESeq2 v1.18.1. Pairwise comparison or one-way analysis of 
#'variance were run with a parametric fit and genotype as the source 
#'of variation (factor: ‘mutant’ or ‘control’)
Genome_b


