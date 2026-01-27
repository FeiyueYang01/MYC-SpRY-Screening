#!/bin/bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/STAR_v2.7.10b

input_dir="/media/protein/fy2306/projects/base_editing/data/ASO_rnaseq/fastp/batch_2"
output_dir="/media/protein/fy2306/projects/base_editing/data/ASO_rnaseq/bam/batch_2"
STAR_INDEX="/media/scratch/fy2306/genomes/hg38/STAR_v2.7.10b"

# samples=("XP047" "XP048" "XP049" "XP050" "XP051" "XP052" "XP053" "XP054" "XP055" "XP056" "XP057" "XP058" "XP059" "XP060")
samples=("XP059" "XP060")
for sample in "${samples[@]}"; do
	STAR --runThreadN 20 \
	--genomeDir "$STAR_INDEX" \
	--sjdbGTFfile /media/scratch/fy2306/genomes/hg38/GENCODE/gencode.v47.annotation.gtf \
	--readFilesIn "$input_dir/${sample}_R1.fastq.gz" "$input_dir/${sample}_R2.fastq.gz" \
	--readFilesCommand zcat \
	--outFileNamePrefix "$output_dir/${sample}_" \
	--outSAMtype BAM SortedByCoordinate \
	--outBAMsortingThreadN 20 \
	--outFilterMultimapNmax 1 \
	
	echo "Finished mapping $sample..."
done