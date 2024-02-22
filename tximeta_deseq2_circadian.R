# Tximeta - DESeq2 Circadian Project 

#######
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
######

library(devtools)
library(rlang)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
# BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
# BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(ashr)
library(readr)

#Load Tximeta and Tximport
# BiocManager::install('tximeta')
library(tximeta)
# BiocManager::install('tximport')
library(tximport)

#Load DESeq2
# BiocManager::install("DESeq2")
library(DESeq2)

#Load rnaseqGene
# BiocManager::install("rnaseqGene")
library(rnaseqGene)


#Section 1: Creating the dataframe for tximeta 

#Creating a sample table 
#First, list all directories containing data
samples <- list.files(path = "../Salmon/quants/bmal1_batch3", full.names = T, pattern = "quant")
samples

#Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")
files

#Set the name of your files without the excess filename
basenametemp <- tools::file_path_sans_ext(basename(samples))
basename <- gsub(pattern = "\\_quant$", "", basenametemp)
basename

#To check that this has worked: expected outcome is True for each file. 
file.exists(files)

#To add the conditions to your samples
cond <- read.csv("./bmal1/bmal1_batch3.csv", header = TRUE)

#To create a CSV table
coldata <- data.frame(files, names = basename, condition = "", stringsAsFactors = FALSE)
coldata

#Add conditions to coldata column
match_conditions <- match(coldata$names, cond$sample)
coldata$condition <- cond$condition [match_conditions]
coldata


#Linking to local GTF files
indexDir <- file.path("../Salmon/m27index_salmon")
fastaFTP <- c("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.transcripts.fa.gz",
              "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/GRCm39.primary_assembly.genome.fa.gz")
gtfPath <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.primary_assembly.annotation.gtf.gz"
suppressPackageStartupMessages(library(tximeta))
makeLinkedTxome(indexDir=indexDir,
                source="GENCODE",
                organism="Mus musculus",
                release="27",
                genome="GRCm39",
                fasta=fastaFTP,
                gtf=gtfPath,
                write=TRUE)


#Running Tximeta - you need to use to above code, for some reason if i just run the Tximeta on its own,
#it wont link correctly 
se <- tximeta(coldata)

#SummarizedExperiment output
#To remove the files column in our se data
suppressPackageStartupMessages(library(SummarizedExperiment))
colData(se)
#here we show that three matrices that were imported
#you should see counts, abundance, and length were imported
assayNames(se)
#Check if tximeta has imported the correct ranges for the transcripts
#This should be equal to the number of valid targets within the meta_info.json file 
rowRanges(se)
#Also check we have the appropriate genome information, which prevents us from making  mistakes
seqinfo(se)

#Next we need to retrieve the transcript database
#This function lets you access the cached database with the following helper function
edb <- retrieveDb(se)
class(edb)
#The database returned by retrieveDb is either a TxDb in the case of GENCODE or RefSeq GTF annotation file
#Or an EnsDb in the case of an Ensemble GTF annotation file. For further use of these two objects
#Consult the genomicsFeatures vignettes and the ensembldb vignette, respectively (Bioconductor packages)

#Add exons per transcript
#The se object created by tximeta has the start, end, and strand information for each transcript
#Here we swap out the transcript GRanges for exons-by-transcript GRangesList 
#this is a list of GRanges where each element of the list gives the exons for a particular transcript. 
se.exons <- addExons(se)
rowRanges(se.exons)[[1]]
#As with the transcript ranges, the exon ranges will be generated once and cahed locally. 

#Easy summarization to gene-level
#Tximeta can make use of the cached TxDB database for the purpose of summarizing
#transcript-level quantifications and bias corrections to gene-level.
gse <- summarizeToGene(se)
rowRanges(gse)
gse$condition


#####
#THE DESEQDATASET
#Building a DESeqDataSet and begin Differential expression analysis
dds <- DESeqDataSet(gse, design = ~ condition)
dds

#This is essential for diagnostic plots - removes genes with a count less than 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#You need to specify what your conditions are and level them appropriately 
dds$condition <- relevel(dds$condition, ref = "PMWT")
levels(dds$condition)
# dds <- DESeqDataSet(gse, design = ~ condition)

#now you run the DESeq analysis 
dds <- DESeq(dds)
resultsNames(dds)

# results:
res_1 <- results(dds, contrast = c("condition", "AMWT", "PMWT"), alpha = 0.05) 
res_1


# # # CHECKING BATCH EFFECT FIRST 
rld <- rlog(dds)

# # since we are using Gencode, we need to switch the gene version id to the gene id in order to use PCA Explorer.
# #Adding gene names
library(biomaRt)

# #Need to use old mart version 102, otherwise data wont match.
ensembl104 <- useEnsembl(biomart = 'genes',
                         dataset = 'mmusculus_gene_ensembl',
                         version = 104)

genemap <- getBM(attributes = c('ensembl_gene_id_version', 'ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id_version',
                 values = dds@rowRanges$gene_id,
                 mart = ensembl104)

id_match <- match(dds@rowRanges$gene_id, genemap$ensembl_gene_id_version)

dds@rowRanges$gene_id_new <- genemap$ensembl_gene_id [id_match]

# PCA EXplorer
# BiocManager::install("pcaExplorer")
library("pcaExplorer")
genenames <- mapIds(org.Mm.eg.db, keys = dds@rowRanges$gene_id_new, column = "SYMBOL", keytype = "ENSEMBL")
annotation <- data.frame(gene_names = genenames, row.names = rownames(dds), stringsAsFactors = FALSE)
head(annotation)

# Relaunch PCAexplorer to compare the datasets
pcaExplorer(dds = dds, dst = rld, annotation = annotation)



#####
#Annotations
#Adding gene names 
library(biomaRt)
#Next we have to split up the rownames of the results object 
res_1$ensembl_version <- sapply(strsplit(rownames(res_1), split="\\+"), "[",1)
head(res_1,4)

#The next piece of code uses the ENSEMBL mart, querying with Ensembl gene id and requesting
#the entrez gene id and HGNC gene symbol
#Also note: with gencode, i get the ensembl gene id version, so i have to update my biomart info to make sure im 
#getting all the additional info. 

#Need to use old mart version 102, otherwise my data wont match.
ensembl104 <- useEnsembl(biomart = 'genes', 
                         dataset = 'mmusculus_gene_ensembl',
                         version = 104)
# Result 1:
genemap1 <- getBM(attributes = c('ensembl_gene_id_version', 'ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position', 'transcription_start_site', 'strand', 'description'), 
                  filters = 'ensembl_gene_id_version', 
                  values = res_1$ensembl_version, 
                  mart = ensembl104)

id_match1 <- match(res_1$ensembl_version, genemap1$ensembl_gene_id_version)

res_1$ensembl_id <- genemap1$ensembl_gene_id [id_match1]
res_1$gene_name <- genemap1$external_gene_name [id_match1]
res_1$chr <- genemap1$chromosome_name [id_match1]
res_1$start_position <- genemap1$start_position [id_match1]
res_1$end_position <- genemap1$end_position [id_match1]
res_1$tss <- genemap1$transcription_start_site [id_match1]
res_1$strand <- genemap1$strand [id_match1]
res_1$gene_des <- genemap1$description [id_match1]
head(res_1,4)



#Exporting results to CSV files
write.csv(as.data.frame(res_1),
          file="condition_AMWTvPMWT_new.csv")
#DONE

#Easy summarization to gene-level
# Export gene level summary
write.csv(as.data.frame(gse@assays@data@listData),
          file="./bmal1/WT_only/gse.csv")
# gene name
write.csv(as.data.frame(gse@rowRanges@elementMetadata@listData$gene_name),
          file="./bmal1/WT_only/gene_name.csv")


#If gene Id doesnt come up:
genemap_gse <- getBM(attributes = c('ensembl_gene_id_version', 'ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id_version',
                 values = gse@rowRanges@elementMetadata@listData$gene_id,
                 mart = ensembl104)

id_match2 <- match(gse@rowRanges@elementMetadata@listData$gene_id, genemap_gse$ensembl_gene_id_version)

gse@rowRanges@elementMetadata@listData$gene_name <- genemap_gse$external_gene_name [id_match2]

# Merging the quant files of the included data (aka removing the same batch effect issues as Matt had)

library(devtools)
library(rlang)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(ashr)
library(readr)


Bmal1=read.csv('./correlation/rna_fpkm.csv')
Lithium=read.csv('./correlation/atac_fpkm.csv')

gfp=read.csv('./correlation/atac_gfp.csv')

# correlation - pearson 
res1 <- cor.test(RNA$rna, ATAC$log2.average., 
                 method = "pearson")
res1

res2 <- cor.test(RNA$rna, ATAC$log2.average., 
                 method = "spearman")
res2


library("ggpubr")
ggscatter(ATAC, x = "log2.average.", y = "rna", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ATAC", ylab = "RNA")

quantile(ATAC$log2.average., probs = c(0.5, 0.70))
quantile(ATAC$rna, probs = c(0.5, 0.70))


ggscatter(ATAC, x = "log2.average.", y = "rna", color = "Groups",
          palette = c("#882255", "#117733", "#C9D5B5", "#332288","#000000"),
          xlab = "ATAC", ylab = "RNA")



