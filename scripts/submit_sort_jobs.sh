#!/bin/bash

# Define base directories
BASE_BAM_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_output"
SORTED_BAM_BASE_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_output_sorted"

# Create sorted BAM output directory if it doesn't exist
mkdir -p "$SORTED_BAM_BASE_DIR"

# Iterate over each OSD directory
for OSD_DIR in "$BASE_BAM_DIR"/OSD-*; do
    if [ -d "$OSD_DIR" ]; then
        OSD_NAME=$(basename "$OSD_DIR")
        echo "Submitting job for: $OSD_NAME"

        sbatch --export=OSD_NAME="$OSD_NAME",BASE_BAM_DIR="$BASE_BAM_DIR",SORTED_BAM_BASE_DIR="$SORTED_BAM_BASE_DIR" sort_and_index_bam.sh
    fi
done

echo "All jobs submitted."
