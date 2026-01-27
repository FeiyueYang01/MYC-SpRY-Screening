#!/bin/bash

# Cas-OFFinder website can only take 10 sequences at a time for a search with NNN PAM motif  

# Since we are searching without PAM motif restriction, it's just an alignment with x mismatches allowed (not considering gaps)
# Use Bowtie instead
# Bowtie is for short reads (50bp or less) and Bowtie2 is for longer reads  

# fasta file generator: /media/scratch/fy2306/projects/base_editing/script/grna_nc_fasta.ipynb

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/bowtie

# Build the index
# bowtie-build \
# /media/scratch/fy2306/genomes/hg38/GENCODE/GRCh38.primary_assembly.genome.fa \
# /media/scratch/fy2306/genomes/hg38/bowtie/GRCh38.primary_assembly.genome \
# --threads 50
num=2

bowtie \
--threads 50 \
-f \
-v $num \
--all \
--no-unal \
--sam \
/media/scratch/fy2306/genomes/hg38/bowtie/GRCh38.primary_assembly.genome \
/media/scratch/fy2306/projects/base_editing/data/bowtie/input/negative_controls.fa \
/media/scratch/fy2306/projects/base_editing/data/bowtie/negative_controls.bowtie_hg38.mismatch$num.sam \
--un /media/scratch/fy2306/projects/base_editing/data/bowtie/negative_controls.bowtie_hg38.mismatch$num.unmapped.fa \
--al /media/scratch/fy2306/projects/base_editing/data/bowtie/negative_controls.bowtie_hg38.mismatch$num.mapped.fa \