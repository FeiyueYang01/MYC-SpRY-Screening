#!/bin/bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/STAR_v2.7.10b

# output from infer_experiment.py from RSeQC:
# using exon bed as reference
# This is PairEnd Data
# batch 1
# Fraction of reads failed to determine: 0.0725
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.0095
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.9180
# batch 2
# precompiled here: https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/
# This is PairEnd Data
# Fraction of reads failed to determine: 0.3286
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.0091
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.6623

# -a annotation file
# -o output file
# -t feature type, exon by default
# -g attribute type, "gene_id" by default
# -s Perform strand-specific read counting: 0 (unstranded), 1 (stranded) and 2 (reversely stranded, first-strand)
# -p fragments (or templates) will be counted instead of reads.
featureCounts -T 40 -s 2 -p \
-a /media/scratch/fy2306/genomes/hg38/GENCODE/gencode.v47.annotation.gtf \
-o /media/scratch/fy2306/projects/base_editing/data/ASO_rnaseq/batch_2/raw_count_matrix.txt \
/media/protein/fy2306/projects/base_editing/data/ASO_rnaseq/bam/batch_2/*bam

# featureCounts -T 40 -s 2 -p \
# -a /media/scratch/fy2306/genomes/hg38/MANE/release_1.4/MANE.GRCh38.v1.4.refseq_genomic.gtf \
# -o /media/scratch/fy2306/projects/base_editing/data/ASO_rnaseq/raw_count_matrix_MANE.txt \
# /media/protein/fy2306/projects/base_editing/data/ASO_rnaseq/bam/*bam