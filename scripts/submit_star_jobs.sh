#!/bin/bash

# Set the base directory containing trimmed datasets
BASE_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/trimmed"

# Iterate over each trimmed dataset directory and submit a job
for DATASET_DIR in "$BASE_DIR"/OSD-*; do
    if [ -d "$DATASET_DIR" ]; then
        DATASET_NAME=$(basename "$DATASET_DIR")
        echo "Submitting STAR alignment job for: $DATASET_NAME"
        sbatch --export=DATASET_NAME="$DATASET_NAME" star_alignment_job.sh
    fi
done
