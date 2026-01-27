#!/bin/bash
source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/MAGeCK

# https://sourceforge.net/p/mageck/wiki
# https://bitbucket.org/liulab/mageck/src/master/

# tests and ranks sgRNAs and genes based on the read count tables provided
# --day0-label DAY0_LABEL: Specify the label for control sample (usually day 0 or plasmid). For every other sample label, the module will treat it as a treatment condition and compare with control sample.
# --gene-lfc-method
# Method to calculate gene log fold changes (LFC) from sgRNA LFCs. Available methods include the median/mean of all sgRNAs (median/mean), or the median/mean sgRNAs that are ranked in front of the alpha cutoff in RRA (alphamedian/alphamean), or the sgRNA that has the second strongest LFC (secondbest). In the alphamedian/alphamean case, the number of sgRNAs correspond to the "goodsgrna" column in the output, and the gene LFC will be set to 0 if no sgRNA is in front of the alpha cutoff. Default median. (new since v0.5.5)
# --control-gene tells MAGeCK to use provided negative control sgRNAs to generate the null distribution when calculating the p values
method="mean"
win="4"

mageck test \
-k /media/scratch/fy2306/projects/base_editing/data/20250416/mageck/test-1-standard-nucleotide/count/MYC_U1.count_normalized.${win}duplicated.repeats.txt \
--day0-label Cas9-HMOI_D0 \
-n /media/scratch/fy2306/projects/base_editing/data/20250416/mageck/test-1-standard-nucleotide/test/${method}/myc_repeats/${win}/MYC_U1 \
--norm-method none \
--gene-test-fdr-threshold 0.01 \
--adjust-method fdr \
--gene-lfc-method ${method} \
--sort-criteria neg \
--control-sgrna /media/scratch/fy2306/projects/base_editing/data/20250416/mageck/test-1-standard-nucleotide/count/control_sgrna.${win}duplicated.txt \