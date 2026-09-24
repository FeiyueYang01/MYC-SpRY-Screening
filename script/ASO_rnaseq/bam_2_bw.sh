#!/bin/bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/deepTools

bamCoverage -b /media/protein/fy2306/projects/base_editing/data/ASO_rnaseq/bam/batch_2/XP059_Aligned.sortedByCoord.out.bam \
--binSize 5 \
--region chr8 \
--normalizeUsing None \
--centerReads \
-o /media/protein/fy2306/projects/base_editing/data/ASO_rnaseq/bam/batch_2/XP059.chr8.bw

bamCoverage -b /media/protein/fy2306/projects/base_editing/data/ASO_rnaseq/bam/batch_2/XP060_Aligned.sortedByCoord.out.bam \
--binSize 5 \
--region chr8 \
--normalizeUsing None \
--centerReads \
-o /media/protein/fy2306/projects/base_editing/data/ASO_rnaseq/bam/batch_2/XP060.chr8.bw