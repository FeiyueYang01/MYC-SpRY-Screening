#!/usr/bin/env bash

unset LD_LIBRARY_PATH
set -euo pipefail

# BAM filenames or full paths
SAMPLES=(
    "D11_ctl"
    "D11_SCC2"
    "D11_SCC5"
    "MB231_ctl"
    "MB231_ASO723"
    "MB231_ASO723-2"
)

REPLICATES=(R1 R2)
LIBRARIES=(IP input)

INPUT_DIR="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721/bigwigs/tmp"
OUTPUT_DIR="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721-bl/bigwigs/tmp"
BLACKLIST="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721/ENCFF356LFX.bed"
THREADS=16

mkdir -p "$OUTPUT_DIR"

for sample in "${SAMPLES[@]}"; do
    for replicate in "${REPLICATES[@]}"; do
        for library in "${LIBRARIES[@]}"; do

            prefix="${sample}_${replicate}.${library}"
            input_bam="${INPUT_DIR}/${prefix}.human.bam"
            output_bam="${OUTPUT_DIR}/${prefix}.human.filteredbl.bam"
            blacklist_names="${OUTPUT_DIR}/${prefix}.blacklist_qnames.txt"

            echo "Processing: $prefix"

            if [[ ! -s "$input_bam" ]]; then
                echo "WARNING: missing or empty BAM: $input_bam" >&2
                continue
            fi

            # Find read names where either mate overlaps a blacklist region
            bedtools intersect \
				-u \
                -abam "$input_bam" \
                -b "$BLACKLIST" |
                samtools view - |
                cut -f1 |
                LC_ALL=C sort -u \
                > "$blacklist_names"

            # Remove all records with those read names
            if [[ -s "$blacklist_names" ]]; then
                samtools view \
                    -@ "$THREADS" \
                    -b \
                    -N "^${blacklist_names}" \
                    "$input_bam" |
                    samtools sort \
                        -@ "$THREADS" \
                        -o "$output_bam" \
                        -
            else
                echo "No blacklist-overlapping reads found."

                samtools view \
                    -@ "$THREADS" \
                    -b \
                    "$input_bam" |
                    samtools sort \
                        -@ "$THREADS" \
                        -o "$output_bam" \
                        -
            fi

            samtools index -@ "$THREADS" "$output_bam"
			before=$(samtools view -c -f 64 "$input_bam")
			after=$(samtools view -c -f 64 "$output_bam")

			awk \
				-v name="$prefix" \
				-v before="$before" \
				-v after="$after" \
				'BEGIN {
					removed = before - after
					percent = before > 0 ? 100 * removed / before : 0

					printf "%s: removed %d of %d fragments (%.2f%%)\n",
						name, removed, before, percent
				}'

            rm -f "$blacklist_names"

            echo "Created: $output_bam"
            echo
        done
    done
done