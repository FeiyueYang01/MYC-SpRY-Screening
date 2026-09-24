#!/bin/bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/STAR_v2.7.10b

input_dir="/media/protein/fy2306/projects/base_editing/data/PMID_35561311/fastp"
output_dir="/media/protein/fy2306/projects/base_editing/data/PMID_35561311/bam"
STAR_INDEX="/media/dna/fy2306/genomes/hg38/STAR_v2.7.10b"

samples=("SRR17933302" "SRR17933303" "SRR17933304" "SRR17933305")
for sample in "${samples[@]}"; do
	STAR --runThreadN 20 \
	--genomeDir "$STAR_INDEX" \
	--sjdbGTFfile /media/dna/fy2306/genomes/hg38/GENCODE/gencode.v47.annotation.gtf \
	--readFilesIn "$input_dir/${sample}_R1.clean.fastq.gz" "$input_dir/${sample}_R2.clean.fastq.gz" \
	--readFilesCommand zcat \
	--outFileNamePrefix "$output_dir/${sample}_" \
	--outSAMtype BAM SortedByCoordinate \
	--outBAMsortingThreadN 20 \
	--outFilterMultimapNmax 1 \
	
	echo "Finished mapping $sample..."
done