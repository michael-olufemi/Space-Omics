#!/bin/bash

# Define paths
BASE_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_output"
STRANDEDNESS_FILE="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/infer_experiment/strandedness_summary.txt"
RSEM_REF="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_genome/RSEM_ref"
OUTPUT_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/RSEM_counts"
THREADS=32  # Adjust as needed

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Loop through each project directory in STAR_output
for PROJECT_DIR in "${BASE_DIR}"/*; do
    if [[ -d "$PROJECT_DIR" ]]; then
        PROJECT_NAME=$(basename "$PROJECT_DIR")
        PROJECT_OUTPUT_DIR="${OUTPUT_DIR}/${PROJECT_NAME}"

        # Create project-specific output directory
        mkdir -p "${PROJECT_OUTPUT_DIR}"

        # Generate batch job script
        JOB_SCRIPT="${PROJECT_OUTPUT_DIR}/run_rsem_${PROJECT_NAME}.sh"
        cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=RSEM_${PROJECT_NAME}
#SBATCH --output=${PROJECT_OUTPUT_DIR}/rsem_%j.out
#SBATCH --error=${PROJECT_OUTPUT_DIR}/rsem_%j.err
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=500G  # Memory allocation
#SBATCH -p cpu  # Partition
#SBATCH --qos=long  # Quality of service
#SBATCH -t 6-00:00:00  # Job runtime (1 day)

conda activate rsem_env
module load conda/latest
echo "Processing project: ${PROJECT_NAME}"

# Loop through each BAM file in the project directory
for BAM_FILE in ${PROJECT_DIR}/*_Aligned.toTranscriptome.out.bam; do
    SAMPLE_NAME=\$(basename "\$BAM_FILE" _Aligned.toTranscriptome.out.bam)

    # Get the strandedness from the summary file
    STRANDEDNESS=\$(grep "\$SAMPLE_NAME" "$STRANDEDNESS_FILE" | awk '{print \$3}')

    if [[ -z "\$STRANDEDNESS" ]]; then
        echo "Warning: No strandedness found for \$SAMPLE_NAME. Defaulting to 'none'."
        STRANDEDNESS="none"
    fi

    echo "Processing \$SAMPLE_NAME with strandedness: \$STRANDEDNESS"

    rsem-calculate-expression --num-threads ${THREADS} \\
        --alignments \\
        --bam \\
        --paired-end \\
        --seed 12345 \\
        --seed-length 20 \\
        --estimate-rspd \\
        --no-bam-output \\
        --strandedness \$STRANDEDNESS \\
        "\$BAM_FILE" \\
        "${RSEM_REF}" \\
        "${PROJECT_OUTPUT_DIR}/\$SAMPLE_NAME"
done

echo "RSEM processing completed for ${PROJECT_NAME}"
EOF

        # Submit the batch job
        sbatch "$JOB_SCRIPT"
        echo "Submitted job for ${PROJECT_NAME}"
    fi
done
