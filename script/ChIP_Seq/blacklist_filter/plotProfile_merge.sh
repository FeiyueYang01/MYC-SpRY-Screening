#!/usr/bin/env bash

source /media/scratch/fy2306/miniconda3/bin/activate /media/scratch/fy2306/miniconda3/envs/deepTools
unset LD_LIBRARY_PATH

INPUT_DIR="/media/protein/fy2306/projects/base_editing/data/MYC-ChIP-20260721-bl/"
SUFFIX="merged.IP.filteredbl.CPM.equal_depth.bw"

computeMatrix reference-point \
-S "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/D11_ctl.${SUFFIX}" "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/D11_SCC2.${SUFFIX}" "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/D11_SCC5.${SUFFIX}" \
-R "${INPUT_DIR}/peaks/D11_ctl_equal_depth/D11_ctl_equal_depth_peaks.narrowPeak" \
--samplesLabel "Ctl" "SCC2" "SCC5" \
-a 2000 \
-b 2000 \
--referencePoint center \
--binSize 50 \
-p 20 \
-o "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/D11.matrix.gz"

plotProfile -m "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/D11.matrix.gz" \
-o "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/D11.profile.pdf" \
--perGroup \
--plotHeight 10 \
--plotWidth 10 \
--regionsLabel "D11" \
--colors "#C90050" "#258FC7" "#107F80"

computeMatrix reference-point \
-S "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/MB231_ctl.${SUFFIX}" "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/MB231_ASO723.${SUFFIX}" "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/MB231_ASO723-2.${SUFFIX}" \
-R "${INPUT_DIR}/peaks/MB231_ctl_equal_depth/MB231_ctl_equal_depth_peaks.narrowPeak" \
--samplesLabel "Ctl" "ASO723" "ASO723-2" \
-a 2000 \
-b 2000 \
--referencePoint center \
--binSize 50 \
-p 20 \
-o "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/MB231.matrix.gz"

plotProfile -m "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/MB231.matrix.gz" \
-o "${INPUT_DIR}/bigwigs_CPM_merged_equal_depth/MB231.profile.pdf" \
--perGroup \
--plotHeight 10 \
--plotWidth 10 \
--regionsLabel "MDA-MB-231" \
--colors "#C90050" "#258FC7" "#107F80"