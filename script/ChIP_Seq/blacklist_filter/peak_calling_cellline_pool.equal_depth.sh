#!/usr/bin/env bash

# https://support.bioconductor.org/p/108458/
# https://support.bioconductor.org/p/105098/

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/MACS3
unset LD_LIBRARY_PATH

set -euo pipefail

INPUT_DIR="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721-bl/bigwigs/tmp"
OUTPUT_DIR="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721-bl/peaks"

PVALUE="0.05"
THREADS="${THREADS:-16}"
SEED=42

D11_IP_BAMS=(
    "${INPUT_DIR}/D11_ctl_R1.IP.human.filteredbl.bam"
    "${INPUT_DIR}/D11_ctl_R2.IP.human.filteredbl.bam"
    "${INPUT_DIR}/D11_SCC2_R1.IP.human.filteredbl.bam"
    "${INPUT_DIR}/D11_SCC2_R2.IP.human.filteredbl.bam"
    "${INPUT_DIR}/D11_SCC5_R1.IP.human.filteredbl.bam"
    "${INPUT_DIR}/D11_SCC5_R2.IP.human.filteredbl.bam"
)

D11_INPUT_BAMS=(
    "${INPUT_DIR}/D11_ctl_R1.input.human.filteredbl.bam"
    "${INPUT_DIR}/D11_ctl_R2.input.human.filteredbl.bam"
    "${INPUT_DIR}/D11_SCC2_R1.input.human.filteredbl.bam"
    "${INPUT_DIR}/D11_SCC2_R2.input.human.filteredbl.bam"
    "${INPUT_DIR}/D11_SCC5_R1.input.human.filteredbl.bam"
    "${INPUT_DIR}/D11_SCC5_R2.input.human.filteredbl.bam"
)

MB231_IP_BAMS=(
    "${INPUT_DIR}/MB231_ctl_R1.IP.human.filteredbl.bam"
    "${INPUT_DIR}/MB231_ctl_R2.IP.human.filteredbl.bam"
    "${INPUT_DIR}/MB231_ASO723_R1.IP.human.filteredbl.bam"
    "${INPUT_DIR}/MB231_ASO723_R2.IP.human.filteredbl.bam"
    "${INPUT_DIR}/MB231_ASO723-2_R1.IP.human.filteredbl.bam"
    "${INPUT_DIR}/MB231_ASO723-2_R2.IP.human.filteredbl.bam"
)

MB231_INPUT_BAMS=(
    "${INPUT_DIR}/MB231_ctl_R1.input.human.filteredbl.bam"
    "${INPUT_DIR}/MB231_ctl_R2.input.human.filteredbl.bam"
    "${INPUT_DIR}/MB231_ASO723_R1.input.human.filteredbl.bam"
    "${INPUT_DIR}/MB231_ASO723_R2.input.human.filteredbl.bam"
    "${INPUT_DIR}/MB231_ASO723-2_R1.input.human.filteredbl.bam"
    "${INPUT_DIR}/MB231_ASO723-2_R2.input.human.filteredbl.bam"
)

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 ||
        die "Required command not found: $1"
}

count_fragments() {
    local bam="$1"

    # Count one alignment per paired-end fragment:
    #   -f 66: read 1 from a proper pair
    #   -F 2304: exclude secondary and supplementary alignments
    samtools view \
        -@ "${THREADS}" \
        -c \
        -f 66 \
        -F 2304 \
        "${bam}"
}

min_depth() {
    local min=""
    local depth
    local bam

    for bam in "$@"; do
        [[ -f "${bam}" ]] || die "BAM file not found: ${bam}"

        depth="$(count_fragments "${bam}")"
        if [[ -z "${min}" || "${depth}" -lt "${min}" ]]; then
            min="${depth}"
        fi
    done

    printf "%s\n" "${min}"
}

downsample_fraction() {
    local target_depth="$1"
    local source_depth="$2"

    awk \
        -v target="${target_depth}" \
        -v source="${source_depth}" \
        'BEGIN {
            if (source <= 0) {
                exit 1
            }
            fraction = target / source
            if (fraction > 1) {
                fraction = 1
            }
            printf "%.9f", fraction
        }'
}

downsample_group() {
    local label="$1"
    local target_depth="$2"
    local downsample_dir="$3"
    shift 3

    local downsampled_bams=()
    local counts_file="${downsample_dir}/${label}.downsample_counts.tsv"
    local bam
    local sample_name
    local output_bam
    local source_depth
    local fraction
    local output_depth
	local sample_index=0
	local sample_seed

    printf "sample\tsource_bam\tdownsampled_bam\ttarget_fragments\tsource_fragments\tfraction\toutput_fragments\n" \
        > "${counts_file}"

    for bam in "$@"; do
        [[ -f "${bam}" ]] || die "BAM file not found: ${bam}"

        sample_name="$(basename "${bam}" .bam)"
        output_bam="${downsample_dir}/${sample_name}.equal_depth.bam"
        source_depth="$(count_fragments "${bam}")"
        fraction="$(downsample_fraction "${target_depth}" "${source_depth}")" ||
            die "Cannot downsample ${bam}; source depth is ${source_depth}"
		sample_seed=$((SEED + sample_index))
		sample_index=$((sample_index + 1))

        echo "Downsampling ${sample_name}: ${source_depth} -> ${target_depth} fragments (${fraction})" >&2

		samtools view \
			-@ "${THREADS}" \
			--subsample-seed "${sample_seed}" \
			--subsample "${fraction}" \
			-b \
			-o "${output_bam}" \
			"${bam}"

        samtools index -@ "${THREADS}" "${output_bam}"
        output_depth="$(count_fragments "${output_bam}")"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${sample_name}" \
            "${bam}" \
            "${output_bam}" \
            "${target_depth}" \
            "${source_depth}" \
            "${fraction}" \
            "${output_depth}" \
            >> "${counts_file}"

        downsampled_bams+=("${output_bam}")
    done

    printf "%s\n" "${downsampled_bams[@]}"
}

call_cellline_peaks() {
    local cellline="$1"
    local peak_name="$2"
    local output_dir="$3"
    shift 3

    local -n ip_bams_ref="$1"
    local -n input_bams_ref="$2"
    local downsample_dir="${output_dir}/bam_downsample"
    local ip_target_depth
    local input_target_depth
    local -a downsampled_ip_bams
    local -a downsampled_input_bams

    mkdir -p "${downsample_dir}"

    echo "Counting ${cellline} IP fragments"
    ip_target_depth="$(min_depth "${ip_bams_ref[@]}")"

    echo "Counting ${cellline} input fragments"
    input_target_depth="$(min_depth "${input_bams_ref[@]}")"

    echo "${cellline} IP target depth: ${ip_target_depth} fragments"
    echo "${cellline} input target depth: ${input_target_depth} fragments"

    mapfile -t downsampled_ip_bams < <(
        downsample_group "${cellline}.IP" "${ip_target_depth}" "${downsample_dir}" "${ip_bams_ref[@]}"
    )
    mapfile -t downsampled_input_bams < <(
        downsample_group "${cellline}.input" "${input_target_depth}" "${downsample_dir}" "${input_bams_ref[@]}"
    )

    echo "Calling pooled equal-depth peaks for ${cellline}"

    macs3 callpeak \
        --treatment "${downsampled_ip_bams[@]}" \
        --control "${downsampled_input_bams[@]}" \
        --format BAMPE \
        --name "${peak_name}" \
        --outdir "${output_dir}" \
        --pvalue "${PVALUE}" \
        2>&1 | tee "${output_dir}/${peak_name}.macs3.log"
}

require_command samtools
require_command macs3
require_command awk

mkdir -p "${OUTPUT_DIR}/D11_pool_equal_depth"
mkdir -p "${OUTPUT_DIR}/MB231_pool_equal_depth"

call_cellline_peaks \
    "D11" \
    "D11_pool_equal_depth" \
    "${OUTPUT_DIR}/D11_pool_equal_depth" \
    D11_IP_BAMS \
    D11_INPUT_BAMS

call_cellline_peaks \
    "MB231" \
    "MB231_pool_equal_depth" \
    "${OUTPUT_DIR}/MB231_pool_equal_depth" \
    MB231_IP_BAMS \
    MB231_INPUT_BAMS
