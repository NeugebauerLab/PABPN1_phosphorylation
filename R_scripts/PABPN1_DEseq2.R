#Script for DESeq2 analysis of PABPN1 mutants from poly(A)+ RNA-seq

# load_packages -----------------------------------------------------------
#Load the DEseq2 and tibble packages
library(DESeq2)
library(tibble)
library(pheatmap)
library(stringr)
library(dplyr)
library(EnhancedVolcano)
library(EnsDb.Hsapiens.v79)
library(tximport)
library(viridis)

# load_dds ----------------------------------------------------------------

###Create dds -> This does not need to be changed for each comparison
#Set working directory for the analysis 
setwd("./")

# locate files from an input table of sample information
dir <- "path/to/salmon/quants/"
samples <- read.table("salmon/samples.txt", header=TRUE)
files <- file.path(dir, samples$file_name)
names(files) <- samples$sample

# create db of gene names from ensembl to match gene id
edb <- EnsDb.Hsapiens.v79
genesdb <- genes(edb, return.type="DataFrame")
gene.names <- data.frame(genesdb$gene_id, genesdb$gene_name, genesdb$gene_biotype)
names(gene.names) <- c("row", "symbol", "type")

# create a table that maps transcript ID to gene name
txdb <- transcripts(edb, return.type="DataFrame")


tx2gene <- data.frame(txdb$tx_id, txdb$gene_id)
colnames(tx2gene) <- c("TXNAME", "GENENAME")

# import transcript abundance files
txi <- tximport(files, type="salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
#write.table(txi$abundance,"salmon_exp_table.txt", quote = FALSE, col.names = TRUE, row.names = TRUE)

# create a DESeq object
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ condition)

# filter rows to keep only transcripts with >= 10 read total
keep <- rowSums(counts(dds)) >= 10
sum(keep)
dds <- dds[keep, ]

#PCA of gene expression for replicates
pcaData <- vst(dds, blind=TRUE) %>%  
            plotPCA(intgroup=c("condition", "Replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#plot PCA, color by condition
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# create_results_tables ---------------------------------------------------

# run DE analysis 
dds <- DESeq(dds)

# Generate results - use contrast to determine the conditions
res_EV_WT = results(dds, tidy=TRUE, contrast=c("condition", "h112", "h113"))
res_EV_WT <- merge(res_EV_WT, gene.names, all.x=TRUE)
write.table(res_EV_WT,"Diff_exp_resdf_EV_WT.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

res_WT_4SA = results(dds, tidy=TRUE, contrast=c("condition", "h113", "h115"))
res_WT_4SA <- merge(res_WT_4SA, gene.names, all.x=TRUE)
write.table(res_WT_4SA,"Diff_exp_resdf_WT_4SA.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

res_WT_4SD = results(dds, tidy=TRUE, contrast=c("condition", "h113", "h117"))
res_WT_4SD <- merge(res_WT_4SD, gene.names, all.x=TRUE)
write.table(res_WT_4SD,"Diff_exp_resdf_WT_4SD.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

res_4SA_4SD = results(dds, tidy=TRUE, contrast=c("condition", "h115", "h117"))
res_4SA_4SD <- merge(res_4SA_4SD, gene.names, all.x=TRUE)
write.table(res_4SA_4SD,"Diff_exp_resdf_4SA_4SD.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)
