#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=200G  # Memory allocation
#SBATCH -p cpu  # Partition
#SBATCH --qos=long  # Quality of service
#SBATCH -t 1-00:00:00  # Job runtime (1 day)
#SBATCH -o infer_experiment_${OSD_NAME}.log  # Output log file (per directory)


# Define directories (passed from sbatch)
BAM_DIR="${SORTED_BAM_BASE_DIR}/${OSD_NAME}"
OUTPUT_DIR="${INFER_EXPT_OUTPUT_DIR}/${OSD_NAME}"

# Create output directory if not exists
mkdir -p "$OUTPUT_DIR"

echo "Processing directory: $BAM_DIR"

# Process each sorted BAM file in the OSD directory
for BAM_FILE in "$BAM_DIR"/*Aligned.sortedByCoord_sorted.out.bam; do
    # Extract sample name
    SAMPLE_NAME=$(basename "$BAM_FILE" | sed 's/_Aligned.sortedByCoord_sorted.out.bam//')

    echo "Running infer_experiment.py for sample: $SAMPLE_NAME"

    # Run infer_experiment.py
    infer_experiment.py -r "$BED_FILE" \
        -i "$BAM_FILE" \
        -s 15000000 > "$OUTPUT_DIR/${SAMPLE_NAME}_infer_expt.out"

    echo "Completed infer_experiment.py for $SAMPLE_NAME"
done

echo "Inference experiment analysis completed for directory: $OSD_NAME"
