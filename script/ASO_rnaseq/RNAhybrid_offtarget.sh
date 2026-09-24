#!/bin/bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/RNAhybrid

# We selected the following settings: maximum internal loop length and bulge size were set to 5 nt each, and the free energy threshold was variable depending on the ASO, but was at least 10 kcal/mol lower than the on-target binding free energy. We did not set helix constraints, and chose the compact output to help with data formatting.
RNAhybrid -c -u 5 -v 5 -s 3utr_human -p 1 \
-t /media/scratch/fy2306/genomes/hg38/GENCODE/GENCODE_v47.exon200nt.fasta \
-q /media/scratch/fy2306/projects/base_editing/data/ASO_rnaseq/ASOs.fa \
> /media/scratch/fy2306/projects/base_editing/data/ASO_rnaseq/ASO_offtarget/RNAhybrid.out

