# initialize
rm(list=ls()) ; gc();
.libPaths('/media/scratch/fy2306/miniconda3/envs/R4.4/lib/R/library')

setwd('/media/scratch/fy2306/projects/base_editing/results/PMID_35561311')
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
raw_count_path <- '/media/scratch/fy2306/projects/base_editing/data/PMID_35561311/raw_count_matrix.txt'
raw_count_data <- read.table(raw_count_path, header = TRUE, sep = "\t", comment.char = "#")
raw_count_data <- raw_count_data[, !(names(raw_count_data) %in% c("Chr", "Start", "End", "Strand", "Length"))]
colnames(raw_count_data) <- c("Geneid", "SRR17933302", "SRR17933303", "SRR17933304", "SRR17933305")
raw_count_data <- raw_count_data %>% remove_rownames %>% column_to_rownames(var="Geneid")
head(raw_count_data)

metadata_path <- '/media/scratch/fy2306/projects/base_editing/data/PMID_35561311/metadata.tsv'
metadata <- read.table(metadata_path, header = TRUE, sep = "\t")
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="sample")
metadata$condition <- factor(metadata$condition, levels = c("STING", "GFP"))
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
  scale_shape_manual(values = 1:4)
ggsave("/media/scratch/fy2306/projects/base_editing/results/PMID_35561311/deseq2_pca.png", plot=p2, width=6, height=6)

# quick check of what we have
as.data.frame(colData(dds))
dds$condition <- relevel(dds$condition, ref="GFP")

# deseq2
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("condition", "STING", "GFP")))
res$gene_id <- row.names(res)
summary(res)

# add gene_name and save
mapping_path <- '/media/dna/fy2306/genomes/hg38/GENCODE/gencode.v47.annotation.gene_id2gene_name.tsv'
mapping <- read.table(mapping_path, sep = "\t", header = FALSE, col.names = c("gene_id", "gene_name"))
res <- merge(res, mapping, by = "gene_id", all.x = TRUE)
res$gene_id <- sub("\\..*", "", res$gene_id)
write.table(res, "/media/scratch/fy2306/projects/base_editing/results/PMID_35561311/de_deseq2.tsv", sep="\t", row.names = FALSE, quote=FALSE)

converted_ids <- bitr(res$gene_id,
                  fromType = "ENSEMBL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

res <- left_join(res, converted_ids, by = c("gene_id" = "ENSEMBL"))
res_clean <- res %>% filter(!is.na(ENTREZID))
# If there are duplicates, choose the one with highest absolute LFC
res_clean <- res_clean %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

foldchanges <- res_clean$stat	# or log2FoldChange
names(foldchanges) <- res_clean$ENTREZID
foldchanges <- sort(foldchanges, decreasing = TRUE)

msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    dplyr::rename(ont = gs_name, gene = entrez_gene)
print("msig_h")
table(msig_h$gene %in% names(foldchanges))
length(msig_h$gene)
head(names(foldchanges))

gsea_msig <- GSEA(geneList = foldchanges,
                  TERM2GENE = msig_h,
                  pvalueCutoff = 1,
				  scoreType = "std",
                  verbose = FALSE)

gsea_df <- as.data.frame(gsea_msig)
write.table(gsea_df, file = "/media/scratch/fy2306/projects/base_editing/results/PMID_35561311/gsea_msig_H_std.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
