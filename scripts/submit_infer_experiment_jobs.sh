#!/bin/bash

# Define base directories
SORTED_BAM_BASE_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_output_sorted"
INFER_EXPT_OUTPUT_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/infer_experiment"
BED_FILE="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/annotation.bed"

# Create output directory if not exists
mkdir -p "$INFER_EXPT_OUTPUT_DIR"

# Iterate over each OSD directory
for OSD_DIR in "$SORTED_BAM_BASE_DIR"/OSD-*; do
    if [ -d "$OSD_DIR" ]; then
        OSD_NAME=$(basename "$OSD_DIR")
        echo "Submitting job for: $OSD_NAME"

        sbatch --export=OSD_NAME="$OSD_NAME",SORTED_BAM_BASE_DIR="$SORTED_BAM_BASE_DIR",INFER_EXPT_OUTPUT_DIR="$INFER_EXPT_OUTPUT_DIR",BED_FILE="$BED_FILE" run_infer_experiment.sh
    fi
done

echo "All jobs submitted."
