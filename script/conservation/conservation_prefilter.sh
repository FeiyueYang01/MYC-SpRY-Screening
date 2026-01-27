#!/bin/bash

# wigFix: fixedStep chrom=chr8 start=60001 step=1

gzip -dc /media/scratch/fy2306/projects/base_editing/data/conservation/chr8.phyloP30way.wigFix.gz | awk '
BEGIN {current_chrom = ""; start = 0; step = 0;}
/^fixedStep/ {
    # Parse the header
    split($0, arr, " ");
    for (i in arr) {
        if (arr[i] ~ /^chrom=/) chrom = substr(arr[i], 7);
        if (arr[i] ~ /^start=/) start = substr(arr[i], 7);
        if (arr[i] ~ /^step=/) step = substr(arr[i], 6);
    }
    current_chrom = chrom;
    current_start = start;
    current_step = step;
}
/^-?[0-9]+(\.[0-9]+)?$/ {
    # Only print scores in the desired region
    if (current_chrom == "chr8" && current_start >= 127734052 && current_start <= 127744108) {
        print current_chrom "\t" current_start "\t" $1;
    }
    current_start += current_step;  # Increment position
}' > /media/scratch/fy2306/projects/base_editing/data/conservation/chr8_myc.phyloP30way.wigFix