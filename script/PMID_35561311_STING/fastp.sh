#!/bin/bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/STAR_v2.7.10b

input_dir="/media/protein/fy2306/projects/base_editing/data/PMID_35561311/fastq"
output_dir="/media/protein/fy2306/projects/base_editing/data/PMID_35561311/fastp"
mkdir -p "$output_dir"

samples=("SRR17933302" "SRR17933303" "SRR17933304" "SRR17933305")

for sample in "${samples[@]}"; do
    r1="$input_dir/${sample}_1.fastq.gz"
    r2="$input_dir/${sample}_2.fastq.gz"

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