#!/usr/bin/env bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/MACS3
unset LD_LIBRARY_PATH

set -euo pipefail

INPUT_DIR="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721-bl-D11/peaks"
OUTPUT_DIR="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721-bl-D11/peaks"

QVALUE="0.05"

SAMPLES=(
    # "MB231_ctl"
	# "MB231_ASO723"
	# "MB231_ASO723-2"
	"D11_ctl"
	"D11_SCC2"
	"D11_SCC5"
)

mkdir -p "${OUTPUT_DIR}"

for sample in "${SAMPLES[@]}"; do

    cellline="${sample%%_*}"
    ip_r1="${INPUT_DIR}/${cellline}_pool_equal_depth/bam_downsample/${sample}_R1.IP.human.filteredbl.equal_depth.bam"
    ip_r2="${INPUT_DIR}/${cellline}_pool_equal_depth/bam_downsample/${sample}_R2.IP.human.filteredbl.equal_depth.bam"
    input_r1="${INPUT_DIR}/${cellline}_pool_equal_depth/bam_downsample/${sample}_R1.input.human.filteredbl.equal_depth.bam"
    input_r2="${INPUT_DIR}/${cellline}_pool_equal_depth/bam_downsample/${sample}_R2.input.human.filteredbl.equal_depth.bam"

    for bam in "${ip_r1}" "${ip_r2}" "${input_r1}" "${input_r2}"; do
        [[ -f "${bam}" ]] || {
            echo "ERROR: BAM file not found: ${bam}" >&2
            exit 1
        }
    done

    sample_name="${sample}_equal_depth"
    sample_dir="${OUTPUT_DIR}/${sample_name}"

    mkdir -p "${sample_dir}"

    echo "Calling pooled peaks for ${sample}"

    macs3 callpeak \
        --treatment "${ip_r1}" "${ip_r2}" \
        --control "${input_r1}" "${input_r2}" \
        --format BAMPE \
        --name "${sample_name}" \
        --outdir "${sample_dir}" \
        --qvalue "${QVALUE}" \
        2>&1 | tee "${sample_dir}/${sample_name}.macs3.log"
done

echo "Pooled peak calling complete: ${OUTPUT_DIR}"