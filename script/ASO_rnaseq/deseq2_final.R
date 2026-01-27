# initialize
rm(list=ls()) ; gc();
.libPaths('/media/scratch/fy2306/miniconda3/envs/R4.4/lib/R/library')

setwd('/media/scratch/fy2306/projects/base_editing/results/ASO_rnaseq_batch_2_outliers')
set.seed(1337)

library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)
library(ggplot2)
library(dplyr)

# load in and clean the data
raw_count_path <- '/media/scratch/fy2306/projects/base_editing/data/ASO_rnaseq/batch_2/raw_count_matrix.txt'
raw_count_data <- read.table(raw_count_path, header = TRUE, sep = "\t", comment.char = "#")
raw_count_data <- raw_count_data[, !(names(raw_count_data) %in% c("Chr", "Start", "End", "Strand", "Length"))]
colnames(raw_count_data) <- c("Geneid", "XP047", "XP048", "XP049", "XP050", "XP051", "XP052", "XP053", "XP054", "XP055", "XP056", "XP057", "XP058")
raw_count_data <- raw_count_data %>% remove_rownames %>% column_to_rownames(var="Geneid")
head(raw_count_data)

metadata_path <- '/media/scratch/fy2306/projects/base_editing/data/ASO_rnaseq/batch_2/metadata.tsv'
metadata <- read.table(metadata_path, header = TRUE, sep = "\t")
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="sample")
metadata$condition <- factor(metadata$condition, levels = c("Control_48", "Treatment_48", "Control_72", "Treatment_72"))
head(metadata)

dds <- DESeqDataSetFromMatrix(countData = raw_count_data, colData = metadata, design = ~ condition)

# PCA plot
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p2 <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=name)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_shape_manual(values = 1:12)
ggsave("/media/scratch/fy2306/projects/base_editing/results/ASO_rnaseq_batch_2_outliers/deseq2_pca.png", plot=p2, width=6, height=6)

# keep 48h
samples_48 <- c("XP047","XP048","XP049","XP050","XP051","XP052")
cts48 <- as.matrix(raw_count_data[, samples_48, drop = FALSE])
coldata48 <- metadata[samples_48, , drop = FALSE]
dds_48 <- DESeqDataSetFromMatrix(countData = cts48, colData = coldata48, design = ~ condition)
# quick check of what we have
as.data.frame(colData(dds_48))
dds_48$condition <- relevel(dds_48$condition, ref="Control_48")

# deseq2
smallestGroupSize <- 3
keep <- rowSums(counts(dds_48) >= 10) >= smallestGroupSize
dds_48 <- dds_48[keep,]

dds_48 <- DESeq(dds_48)
res <- as.data.frame(results(dds_48, contrast = c("condition", "Treatment_48", "Control_48")))
res$gene_id <- row.names(res)
summary(res)

# add gene_name and save
mapping_path <- '/media/scratch/fy2306/genomes/hg38/GENCODE/gencode.v47.annotation.gene_id2gene_name.tsv'
mapping <- read.table(mapping_path, sep = "\t", header = FALSE, col.names = c("gene_id", "gene_name"))
res <- merge(res, mapping, by = "gene_id", all.x = TRUE)
res$gene_id <- sub("\\..*", "", res$gene_id)
write.table(res, "/media/scratch/fy2306/projects/base_editing/results/ASO_rnaseq_batch_2_outliers/de_deseq2.48h.tsv", sep="\t", row.names = FALSE, quote=FALSE)


# keep 72h
samples_72 <- c("XP053", "XP054", "XP055", "XP056", "XP057", "XP058")
cts72 <- as.matrix(raw_count_data[, samples_72, drop = FALSE])
coldata72 <- metadata[samples_72, , drop = FALSE]
dds_72 <- DESeqDataSetFromMatrix(countData = cts72, colData = coldata72, design = ~ condition)
# quick check of what we have
as.data.frame(colData(dds_72))
dds_72$condition <- relevel(dds_72$condition, ref="Control_72")

# deseq2
smallestGroupSize <- 3
keep <- rowSums(counts(dds_72) >= 10) >= smallestGroupSize
dds_72 <- dds_72[keep,]

dds_72 <- DESeq(dds_72)
res <- as.data.frame(results(dds_72, contrast = c("condition", "Treatment_72", "Control_72")))
res$gene_id <- row.names(res)
summary(res)

# add gene_name and save
mapping_path <- '/media/scratch/fy2306/genomes/hg38/GENCODE/gencode.v47.annotation.gene_id2gene_name.tsv'
mapping <- read.table(mapping_path, sep = "\t", header = FALSE, col.names = c("gene_id", "gene_name"))
res <- merge(res, mapping, by = "gene_id", all.x = TRUE)
res$gene_id <- sub("\\..*", "", res$gene_id)
write.table(res, "/media/scratch/fy2306/projects/base_editing/results/ASO_rnaseq_batch_2_outliers/de_deseq2.72h.tsv", sep="\t", row.names = FALSE, quote=FALSE)