# initialize
rm(list=ls()) ; gc();
.libPaths('/media/scratch/fy2306/miniconda3/envs/R4.4/lib/R/library')

args <- commandArgs(trailingOnly = TRUE)
wdir <- args[1]
setwd(wdir)
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
hour <- args[2]
res_path <- paste0(wdir, '/de_deseq2.', hour, 'h.tsv')
res <- read.table(res_path, header = TRUE, sep = "\t", comment.char = "#")
res_sig <- subset(res, padj < 0.05)	# NAs will be excluded

# GSEA
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

# C2
msig_C2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    dplyr::rename(ont = gs_name, gene = entrez_gene)
print("msig_C2")
table(msig_C2$gene %in% names(foldchanges))

gsea_msig <- GSEA(geneList = foldchanges,
                  TERM2GENE = msig_C2,
                  pvalueCutoff = 1,
				  scoreType = "std",
                  verbose = FALSE)
gsea_df <- as.data.frame(gsea_msig)
set_name <- "DANG_REGULATED_BY_MYC_DN"
p <- gseaplot2(gsea_msig, set_name, title = set_name, base_size=10,  subplots = 1:2, rel_heights = c(1, 0.25))
ggsave(paste0(wdir, "/enrichment_", hour, "h/gsea_msig_C2_std_", set_name, ".pdf"), plot = p, width = 4, height = 3, dpi = 300)
set_name <- "DANG_REGULATED_BY_MYC_UP"
p <- gseaplot2(gsea_msig, set_name, title = set_name, base_size=10,  subplots = 1:2, rel_heights = c(1, 0.25))
ggsave(paste0(wdir, "/enrichment_", hour, "h/gsea_msig_C2_std_", set_name, ".pdf"), plot = p, width = 4, height = 3, dpi = 300)