#!/bin/bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/deepTools

bamCoverage -b /path/to/XP060_Aligned.sortedByCoord.out.bam \
--binSize 5 \
--region chr8 \
--normalizeUsing None \
--centerReads \
-o XP060.chr8.bw