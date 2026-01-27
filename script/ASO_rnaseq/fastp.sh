#!/bin/bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/STAR_v2.7.10b

input_dir="/media/protein/sequencing_data2/250813_XUEBING_PEIGUO_14_HUMAN_RNA_STRDPOLYA_20M_PE75_AVITI/fastq"
output_dir="/media/protein/fy2306/projects/base_editing/data/ASO_rnaseq/fastp/batch_2"
mkdir -p "$output_dir"

samples=("XP047" "XP048" "XP049" "XP050" "XP051" "XP052" "XP053" "XP054" "XP055" "XP056" "XP057" "XP058" "XP059" "XP060")

for sample in "${samples[@]}"; do
    r1="$input_dir/${sample}_R1.fastq.gz"
    r2="$input_dir/${sample}_R2.fastq.gz"

    out_r1="$output_dir/${sample}_R1.clean.fastq.gz"
    out_r2="$output_dir/${sample}_R2.clean.fastq.gz"
    html="$output_dir/${sample}_fastp.html"

	# default: quality filter, length filter, adapter trimming
	# The most widely used adapters are Illumina TruSeq adapters. If your data is from the TruSeq library, fastp should be able to detect it successfully
    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$out_r1" \
        -O "$out_r2" \
        -h "$html" \
        --thread 40

	echo "Finished fastp for $sample..."
done