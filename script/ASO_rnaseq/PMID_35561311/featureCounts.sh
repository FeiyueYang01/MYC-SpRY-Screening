#!/bin/bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/STAR_v2.7.10b

# output from infer_experiment.py from RSeQC:
# using exon bed as reference
# This is PairEnd Data
# Fraction of reads failed to determine: 0.0900
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.4549
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.4551
# unstranded data

# -a annotation file
# -o output file
# -t feature type, exon by default
# -g attribute type, "gene_id" by default
# -s Perform strand-specific read counting: 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
# -p fragments (or templates) will be counted instead of reads.
featureCounts -T 40 -s 0 -p \
-a /media/dna/fy2306/genomes/hg38/GENCODE/gencode.v47.annotation.gtf \
-o /media/scratch/fy2306/projects/base_editing/data/PMID_35561311/raw_count_matrix.txt \
/media/protein/fy2306/projects/base_editing/data/PMID_35561311/bam/*bam