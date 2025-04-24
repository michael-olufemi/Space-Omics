#!/bin/bash

# Set the base directory containing your datasets
BASE_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/downloads2"

# Iterate over each dataset directory and submit a job
for DATASET_DIR in "$BASE_DIR"/OSD-*; do
    if [ -d "$DATASET_DIR" ]; then
        DATASET_NAME=$(basename "$DATASET_DIR")
        sbatch --export=DATASET_NAME="$DATASET_NAME" trim_galore_job.sh
    fi
done
