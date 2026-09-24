#!/usr/bin/env bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/bedtools
unset LD_LIBRARY_PATH

set -euo pipefail

# Convert replicate-level paired-end IP ChIP-seq BAMs to positionally
# deduplicated fragment BED files, then count fragments overlapping pooled
# D11 and MB231 MYC peaks by >5 bp.

DATA_DIR="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721-bl"
BAM_DIR="${DATA_DIR}/bigwigs/tmp"
PEAK_DIR="${DATA_DIR}/peaks"
FRAGMENT_DIR="${DATA_DIR}/fragment_bed"
COUNT_DIR="${DATA_DIR}/peak_fragment_counts"
TMP_BAM_DIR="${FRAGMENT_DIR}/tmp_name_sorted_bams"

MIN_OVERLAP_BP="${MIN_OVERLAP_BP:-5}"
THREADS="${THREADS:-16}"

D11_PEAK_FILE="${PEAK_DIR}/D11_pool_equal_depth/D11_pool_equal_depth_peaks.narrowPeak"
MB231_PEAK_FILE="${PEAK_DIR}/MB231_pool_equal_depth/MB231_pool_equal_depth_peaks.narrowPeak"

D11_COUNT_TSV="${COUNT_DIR}/D11_pool_equal_depth_peaks_fragment_counts.tsv"
MB231_COUNT_TSV="${COUNT_DIR}/MB231_pool_equal_depth_peaks_fragment_counts.tsv"

D11_TMP_COUNT_DIR="${COUNT_DIR}/tmp_D11_pool_counts"
MB231_TMP_COUNT_DIR="${COUNT_DIR}/tmp_MB231_pool_counts"

D11_IP_SAMPLES=(
    "D11_ctl_R1"
    "D11_ctl_R2"
    "D11_SCC2_R1"
    "D11_SCC2_R2"
    "D11_SCC5_R1"
    "D11_SCC5_R2"
)

MB231_IP_SAMPLES=(
    "MB231_ctl_R1"
    "MB231_ctl_R2"
    "MB231_ASO723_R1"
    "MB231_ASO723_R2"
    "MB231_ASO723-2_R1"
    "MB231_ASO723-2_R2"
)

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 ||
        die "Required command not found: $1"
}

join_by_comma() {
    local IFS=,
    echo "$*"
}

make_fragment_bed() {
    local sample="$1"
    local bam="${BAM_DIR}/${sample}.IP.human.filteredbl.bam"
    local name_sorted_bam="${TMP_BAM_DIR}/${sample}.name_sorted.bam"
    local output_bed="${FRAGMENT_DIR}/${sample}.fragments.dedup.bed"
    local tmp_bed="${output_bed}.tmp"

    [[ -f "${bam}" ]] || die "BAM not found for ${sample}: ${bam}"

    echo "Creating positionally deduplicated paired-end fragments for ${sample}"

    samtools sort \
        -@ "${THREADS}" \
        -n \
        -T "${TMP_BAM_DIR}/${sample}.sort" \
        -o "${name_sorted_bam}" \
        "${bam}"

    bedtools bamtobed -bedpe -i "${name_sorted_bam}" |
    awk '
        BEGIN { OFS = "\t" }
        $1 == $4 && $1 != "." {
            start = ($2 < $5) ? $2 : $5
            end = ($3 > $6) ? $3 : $6

            if (end > start) {
                print $1, start, end
            }
        }
    ' |
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n -u > "${tmp_bed}"

    mv "${tmp_bed}" "${output_bed}"
    rm -f "${name_sorted_bam}"
}

count_sample_peaks() {
    local sample="$1"
    local peak_file="$2"
    local tmp_count_dir="$3"
    local fragment_bed="${FRAGMENT_DIR}/${sample}.fragments.dedup.bed"
    local tmp_counts="${tmp_count_dir}/${sample}.counts.tsv"

    [[ -f "${peak_file}" ]] || die "Peak file not found: ${peak_file}"
    [[ -f "${fragment_bed}" ]] || die "Fragment BED not found for ${sample}: ${fragment_bed}"

    echo "Counting ${sample} fragments in $(basename "${peak_file}")"

    bedtools intersect \
        -a "${peak_file}" \
        -b "${fragment_bed}" \
        -wo |
    awk -v min_overlap="${MIN_OVERLAP_BP}" '
        BEGIN { FS = OFS = "\t" }
        $NF > min_overlap {
            key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
            counts[key]++
        }
        END {
            for (key in counts) {
                split(key, fields, SUBSEP)
                print fields[1], fields[2], fields[3], fields[4], counts[key]
            }
        }
    ' > "${tmp_counts}"

    awk \
        'BEGIN { FS = OFS = "\t" }
        NR == FNR {
            key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
            counts[key] = $5
            next
        }
        {
            key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
            count = (key in counts) ? counts[key] : 0
            print $1, $2, $3, $4, count
        }' \
        "${tmp_counts}" \
        "${peak_file}" \
        > "${tmp_counts}.complete"

    mv "${tmp_counts}.complete" "${tmp_counts}"
}

write_count_matrix() {
    local peak_file="$1"
    local tmp_count_dir="$2"
    local output_tsv="$3"
    shift 3
    local samples=("$@")
    local sample_csv
    local count_files=()
    local sample

    sample_csv="$(join_by_comma "${samples[@]}")"

    for sample in "${samples[@]}"; do
        count_files+=("${tmp_count_dir}/${sample}.counts.tsv")
    done

    awk \
        -v peak_file="${peak_file}" \
        -v samples_csv="${sample_csv}" \
        'BEGIN {
            FS = OFS = "\t"
            n_samples = split(samples_csv, samples, ",")
            printf "chrom\tstart\tend\tpeak_id"
            for (i = 1; i <= n_samples; i++) {
                printf "\t%s", samples[i]
            }
            printf "\n"
        }
        FILENAME != peak_file {
            sample = FILENAME
            sub(/^.*\//, "", sample)
            sub(/\.counts\.tsv$/, "", sample)
            key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
            counts[key, sample] = $5
            next
        }
        {
            key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
            printf "%s\t%s\t%s\t%s", $1, $2, $3, $4
            for (i = 1; i <= n_samples; i++) {
                sample = samples[i]
                count = ((key, sample) in counts) ? counts[key, sample] : 0
                printf "\t%s", count
            }
            printf "\n"
        }' \
        "${count_files[@]}" \
        "${peak_file}" \
        > "${output_tsv}"
}

process_peak_set() {
    local peak_file="$1"
    local tmp_count_dir="$2"
    local output_tsv="$3"
    shift 3
    local samples=("$@")
    local sample

    mkdir -p "${tmp_count_dir}"

    for sample in "${samples[@]}"; do
        make_fragment_bed "${sample}"
        count_sample_peaks "${sample}" "${peak_file}" "${tmp_count_dir}"
    done

    write_count_matrix "${peak_file}" "${tmp_count_dir}" "${output_tsv}" "${samples[@]}"
    rm -rf "${tmp_count_dir}"
}

require_command samtools
require_command bedtools
require_command awk
require_command sort

[[ -d "${BAM_DIR}" ]] || die "BAM directory not found: ${BAM_DIR}"
[[ -d "${PEAK_DIR}" ]] || die "Peak directory not found: ${PEAK_DIR}"
[[ -f "${D11_PEAK_FILE}" ]] || die "D11 pooled peak file not found: ${D11_PEAK_FILE}"
[[ -f "${MB231_PEAK_FILE}" ]] || die "MB231 pooled peak file not found: ${MB231_PEAK_FILE}"

mkdir -p "${FRAGMENT_DIR}" "${COUNT_DIR}" "${TMP_BAM_DIR}"

process_peak_set "${D11_PEAK_FILE}" "${D11_TMP_COUNT_DIR}" "${D11_COUNT_TSV}" "${D11_IP_SAMPLES[@]}"
process_peak_set "${MB231_PEAK_FILE}" "${MB231_TMP_COUNT_DIR}" "${MB231_COUNT_TSV}" "${MB231_IP_SAMPLES[@]}"

rm -rf "${TMP_BAM_DIR}"

echo "Wrote paired-end fragment BED files to ${FRAGMENT_DIR}"
echo "Wrote D11 pooled peak count matrix: ${D11_COUNT_TSV}"
echo "Wrote MB231 pooled peak count matrix: ${MB231_COUNT_TSV}"
