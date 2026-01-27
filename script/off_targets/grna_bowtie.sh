#!/bin/bash

# Cas-OFFinder website can only take 10 sequences at a time for a search with NNN PAM motif  
# Cas-OFFinder has a downloadable version as well

# Since we are searching without PAM motif restriction, it's just an alignment with x mismatches allowed (not considering gaps)
# Use Bowtie instead
# Bowtie is for short reads (50bp or less) and Bowtie2 is for longer reads  

# fasta file generator: /media/scratch/fy2306/projects/base_editing/script/off_targets/grna_fasta.ipynb

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/bowtie

# Build the index
# bowtie-build \
# /media/scratch/fy2306/genomes/hg38/GENCODE/GRCh38.primary_assembly.genome.fa \
# /media/scratch/fy2306/genomes/hg38/bowtie/GRCh38.primary_assembly.genome \
# --threads 50

# In the -n alignment mode, an alignment’s “stratum” is defined as the number of mismatches in the “seed” region, i.e. the leftmost L bases, where L is set with the -l option. 
# In the -v alignment mode, an alignment’s stratum is defined as the total number of mismatches in the entire alignment.
# -v 	alignments may have no more than V mismatches, end-to-end alignment
# -f	The query input files (specified either as <m1> and <m2>, or as <s>) are FASTA files
# --all	Report all alignments, not just the best oneReport all valid alignments per read or pair
# --no-unal	Suppress SAM records for reads that failed to align
num=1

bowtie \
--threads 80 \
-f \
-v $num \
--all \
--no-unal \
--sam \
/media/scratch/fy2306/genomes/hg38/bowtie/GRCh38.primary_assembly.genome \
/media/scratch/fy2306/projects/base_editing/data/bowtie/all/input/all_sgrna_seqs.fa \
/media/scratch/fy2306/projects/base_editing/data/bowtie/all/all_sgrna_seqs.bowtie_hg38.mismatch$num.sam \
--un /media/scratch/fy2306/projects/base_editing/data/bowtie/all/all_sgrna_seqs.bowtie_hg38.mismatch$num.unmapped.fa \
--al /media/scratch/fy2306/projects/base_editing/data/bowtie/all/all_sgrna_seqs.bowtie_hg38.mismatch$num.mapped.fa \