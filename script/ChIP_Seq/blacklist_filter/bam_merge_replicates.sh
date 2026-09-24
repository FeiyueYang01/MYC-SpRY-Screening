#!/usr/bin/env bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/samtools
unset LD_LIBRARY_PATH

set -euo pipefail

INPUT_DIR="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721-bl/peaks/MB231_pool_equal_depth/bam_downsample"
OUTPUT_DIR="${INPUT_DIR}/merged"

THREADS=16

SAMPLES=(
    "MB231_ctl"
    "MB231_ASO723"
    "MB231_ASO723-2"
    # "D11_ctl"
    # "D11_SCC2"
    # "D11_SCC5"
)

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 ||
        die "Required command not found: $1"
}

merge_replicates() {
    local sample="$1"
    local assay="$2"
    local r1_bam="${INPUT_DIR}/${sample}_R1.${assay}.human.filteredbl.equal_depth.bam"
    local r2_bam="${INPUT_DIR}/${sample}_R2.${assay}.human.filteredbl.equal_depth.bam"
    local merged_bam="${OUTPUT_DIR}/${sample}.merged.${assay}.human.filteredbl.equal_depth.bam"

    [[ -f "${r1_bam}" ]] || die "Replicate 1 BAM not found: ${r1_bam}"
    [[ -f "${r2_bam}" ]] || die "Replicate 2 BAM not found: ${r2_bam}"

    echo "Merging ${sample} ${assay}"

    samtools merge \
        -@ "${THREADS}" \
        -f \
        "${merged_bam}" \
        "${r1_bam}" \
        "${r2_bam}"

    samtools index -@ "${THREADS}" "${merged_bam}"
}

require_command samtools

mkdir -p "${OUTPUT_DIR}"

for sample in "${SAMPLES[@]}"; do
    merge_replicates "${sample}" "input"
done

echo "Merged BAMs written to ${OUTPUT_DIR}"
