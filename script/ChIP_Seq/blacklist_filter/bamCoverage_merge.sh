#!/usr/bin/env bash

set -euo pipefail

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/deepTools
unset LD_LIBRARY_PATH

# Re-run only the bamCoverage step from bamCoverage.sh.
# This script reuses the human-only BAMs produced by the original script, while
# allowing a different deepTools normalization mode.

BASE_DIR="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721-bl"

NORMALIZE_USING="CPM"
OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/bigwigs_${NORMALIZE_USING}_merged_equal_depth}"

# sample identifiers
SAMPLES=(
    "MB231_ctl"
    "MB231_ASO723"
    "MB231_ASO723-2"
    "D11_ctl"
    "D11_SCC2"
    "D11_SCC5"
)

BIN_SIZE=50
THREADS=16

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 ||
        die "Required command not found: $1"
}

require_bam_index() {
    local bam="$1"

    [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" ]] ||
        die "BAM index not found for reused BAM: ${bam}"
}

require_command bamCoverage

[[ ${#SAMPLES[@]} -gt 0 ]] || die "The SAMPLES list is empty."

mkdir -p "$OUTPUT_DIR"

for sample in "${SAMPLES[@]}"; do
	cellline="${sample%%_*}"
    human_ip_bam="${BASE_DIR}/peaks/${cellline}_pool_equal_depth/bam_downsample/merged/${sample}.merged.IP.human.filteredbl.equal_depth.bam"
    output_bw="${OUTPUT_DIR}/${sample}.merged.IP.filteredbl.${NORMALIZE_USING}.equal_depth.bw"

    [[ -f "$human_ip_bam" ]] ||
        die "Reused human-only IP BAM not found for ${sample}: ${human_ip_bam}"
    require_bam_index "$human_ip_bam"

    echo "Processing ${sample}"
    echo "  Reusing BAM:       ${human_ip_bam}"
    echo "  normalizeUsing:    ${NORMALIZE_USING}"

    bamcoverage_args=(
        --bam "$human_ip_bam"
        --outFileName "$output_bw"
        --outFileFormat bigwig
        --normalizeUsing "$NORMALIZE_USING"
		# --effectiveGenomeSize 2747877702 # Read length 75, GRCh38
        --extendReads
        --samFlagInclude 66
        --binSize "$BIN_SIZE"
        --numberOfProcessors "$THREADS"
    )

    bamCoverage "${bamcoverage_args[@]}"

    echo "  Created: ${output_bw}"
    echo
done

echo "Completed all samples."
echo "Output directory: ${OUTPUT_DIR}"
