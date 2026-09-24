#!/bin/bash
source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/MAGeCK

# https://sourceforge.net/p/mageck/wiki
# https://bitbucket.org/liulab/mageck/src/master/

# /home/xw2629/MYC-screen/ABE-CBE-Cas9-20240430/mageck.sh
# /media/scratch/fy2306/projects/base_editing/script/20241220/mageck_ABEHMOID0.sh

# normalization method:
# https://www.cell.com/cell/fulltext/S0092-8674(21)00084-2
# Read counts for the two sublibraries were independently normalized using the median ratio method. 
# In order to avoid false positives, sgRNAs with low sequencing coverage (average normalized read count at T0 < 200) were eliminated from the analyses. 
# Comparisons of T0 versus T18 were performed using paired MAGeCK robust rank aggregation (RRA) using normalized values as input 
# (e.g., mageck test –k Countfile.txt, -t T18_rep1, T18_rep2, T18_rep3, -c T0_rep1, T0_rep2, T0_rep3 –paired –norm-method none -n T18 –adjust-method fdr).
# Summary: normalized in count, no normalization in test

filelist="/media/protein/sequencing_data2/20250416_AV100007_2025-04-16-BE-SPRY-2ND-Screening/fastq/Cas9-HMOI_D0_R1.fastq.gz \
/media/protein/sequencing_data2/20250416_AV100007_2025-04-16-BE-SPRY-2ND-Screening/fastq/Cas9-LMOI_D20_R1.fastq.gz \
/media/protein/sequencing_data2/20250416_AV100007_2025-04-16-BE-SPRY-2ND-Screening/fastq/Cas9-LMOI_D8_R1.fastq.gz \
/media/rna/sequencing_data/peiguo_shi/20240425_AV100007_2024-04-25-ABE-CBE-SPRY-Screening/fastq/Cas9_SpRY_L_MOI_D_20_R1.fastq.gz \
/media/rna/sequencing_data/peiguo_shi/20240425_AV100007_2024-04-25-ABE-CBE-SPRY-Screening/fastq/Cas9_SpRY_L_MOI_D_8_R1.fastq.gz"

samples=$(for f in $filelist; do basename "$f" | sed 's/_R1\.fastq\.gz$//'; done | paste -sd, -)

# test-1: normalization: count - control
# raw count always stays the same 
# for 2-batch merged samples and separate batch samples: normalized count can be different, but lfc is roughly the same
# using all sgRNAs and MYC only sgRNAs: normalized count stays the same

# collects sgRNA read count information from fastq files
# note that --reverse-complement the library
mageck count \
-l /media/protein/sequencing_data2/20250416_AV100007_2025-04-16-BE-SPRY-2ND-Screening/fastq/MYC-lib-for-mageck.txt \
-n /media/scratch/fy2306/projects/base_editing/data/20250416/mageck/test-1-standard/count/MYC_U1 \
--sample-label $samples \
--fastq $filelist \
--reverse-complement \
--norm-method control \
--control-gene /media/scratch/fy2306/projects/base_editing/data/20240430/mageck/control_gene.txt \

# tests and ranks sgRNAs and genes based on the read count tables provided
# --day0-label DAY0_LABEL: Specify the label for control sample (usually day 0 or plasmid). For every other sample label, the module will treat it as a treatment condition and compare with control sample.
# --control-gene tells MAGeCK to use provided negative control sgRNAs to generate the null distribution when calculating the p values
mageck test \
-k /media/scratch/fy2306/projects/base_editing/data/20250416/mageck/test-1-standard/count/MYC_U1.count_normalized.txt \
--day0-label Cas9-HMOI_D0 \
-n /media/scratch/fy2306/projects/base_editing/data/20250416/mageck/test-1-standard/test/MYC_U1 \
--norm-method none \
--gene-test-fdr-threshold 0.25 \
--adjust-method fdr \
--control-gene /media/scratch/fy2306/projects/base_editing/data/20240430/mageck/control_gene.txt \

# test-2: normalization: count - median
mageck count \
-l /media/protein/sequencing_data2/20250416_AV100007_2025-04-16-BE-SPRY-2ND-Screening/fastq/MYC-lib-for-mageck.txt \
-n /media/scratch/fy2306/projects/base_editing/data/20250416/mageck/test-2-standard/count/MYC_U1 \
--sample-label $samples \
--fastq $filelist \
--reverse-complement \
--norm-method median \

# tests and ranks sgRNAs and genes based on the read count tables provided
# --day0-label DAY0_LABEL: Specify the label for control sample (usually day 0 or plasmid). For every other sample label, the module will treat it as a treatment condition and compare with control sample.
# --control-gene tells MAGeCK to use provided negative control sgRNAs to generate the null distribution when calculating the p values
mageck test \
-k /media/scratch/fy2306/projects/base_editing/data/20250416/mageck/test-2-standard/count/MYC_U1.count_normalized.txt \
--day0-label Cas9-HMOI_D0 \
-n /media/scratch/fy2306/projects/base_editing/data/20250416/mageck/test-2-standard/test/MYC_U1 \
--norm-method none \
--gene-test-fdr-threshold 0.25 \
--adjust-method fdr \
--control-gene /media/scratch/fy2306/projects/base_editing/data/20240430/mageck/control_gene.txt \