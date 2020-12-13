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

Veh_fastq.files <- list.files(path = "rnaseq/", pattern = "^V", full.names = TRUE)
Veh_fastq.files

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
align.stat <- align(index = "Ebi_Homo_Sapien_Index", readfile1 = Veh_fastq.files)



#get annotation file to count mapped reads against
# ann <- getInBuiltAnnotation(annotation = "hg38")

#bild a vector of mapped reads
bam.files <- list.files(path = "rnaseq/", pattern = ".BAM$", full.names = TRUE)


#featureCounts generates a matrix of read counts to genes in each sample
#strandSpecific 'TruSeq Stranded mRNAseq Sample Prep kit' (Illumina). 
#Read 1 aligns to the ANTISENSE strand and Read 2 aligns to the SENSE strand
fc <- featureCounts(bam.files, annot.inbuilt="hg38", strandSpecific = 2)


#This outputs the mapping status of the reads
#fc$stat

#SummarizedExperiment 
#matrix features are accessed with an R function of the same name: assay (counts), rowRanges(geneID) and colData(samples/targets).
se <- SummarizedExperiment(fc$counts, rowRanges = fc$annotation$GeneID, colData = fc$targets)

####write to file
#output a txt file of count data
write.table(x=data.frame(fc$annotation[,c("GeneID")],fc$counts,stringsAsFactors=FALSE),file="strandedcounts.txt",quote=FALSE,sep="\t",row.names=FALSE)

write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],fc$counts,stringsAsFactors=FALSE),file="EbiRscounts.txt",quote=FALSE,sep="\t",row.names=FALSE)
#dimensions - rows are genes, columns are samples
dim(fc)

#load data into a table
EBI_RStrandedCountTable <- read.table(file = 'EbiRscounts.txt', header = TRUE)

#extract null counts as percentage
prop.null <- apply(EBI_RStrandedCountTable, 2, function(x) 100*mean(x==0))
View(prop.null)
print(head(prop.null))
#plot null counts - alot of zeros
barplot(prop.null, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, ylab='Samples', xlab='% of null counts')


#Determine which genes have sufficiently large counts to be retained in a statistical analysis.
#default values of filtByExpr have been used for low expression - 
count.table <- filterByExpr(fc$counts, group = fastq.files)
fc <- fc$counts[count.table,]





#######HISTOGRAM ATTEMPTS ############### FUCK YEAHHHH
ggplot(EBI_RStrandedCountTable) +
  geom_histogram(aes(x = F1_CAAGCTAG.ACATAGCG_L002_R1_001..1..fastq.subread.BAM), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")



hist(counttble$F1_CAAGCTAG.ACATAGCG_L002_R1_001..1..fastq.subread.BAM)




########results
#pca/heatmap/dge list/NB with benjamini-hochberg correction/gene ontology/pattern pathways


