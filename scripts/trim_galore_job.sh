#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=1000G
#SBATCH -p cpu # Partition
#SBATCH --qos=long
#SBATCH -t 6-23:00:00  # Job time limit
#SBATCH -o /scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/logs/trimming_%x.log
#SBATCH --job-name=trim_galore_%x

# Load necessary modules
module load trimgalore/0.6.10
module load cutadapt/3.5-GCCcore-11.2.0

# Set the base directory containing your datasets
BASE_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/downloads2"

# Retrieve the dataset name from the environment variable
DATASET_NAME=${DATASET_NAME}

# Define the dataset directory
DATASET_DIR="${BASE_DIR}/${DATASET_NAME}"

# Define the output directory for trimmed data
OUTPUT_DIR="${BASE_DIR}/trimmed/${DATASET_NAME}"
mkdir -p "$OUTPUT_DIR"

# Change to the dataset directory
cd "$DATASET_DIR" || exit

# Find all FASTQ files in the current dataset directory
FASTQ_FILES=(*_R1_raw.fastq.gz)
NUM_FILES=${#FASTQ_FILES[@]}

if [ $NUM_FILES -eq 0 ]; then
    echo "No FASTQ files found in $DATASET_DIR. Exiting."
    exit 1
fi

# Process each paired-end sample
for R1_FILE in "${FASTQ_FILES[@]}"; do
    # Derive the corresponding R2 filename
    R2_FILE="${R1_FILE/_R1_raw.fastq.gz/_R2_raw.fastq.gz}"

    # Check if the R2 file exists
    if [ ! -f "$R2_FILE" ]; then
        echo "Missing R2 file for $R1_FILE. Skipping this pair."
        continue
    fi

    # Extract sample name (e.g., remove _R1_raw.fastq.gz)
    SAMPLE_NAME="${R1_FILE%_R1_raw.fastq.gz}"
    echo "Trimming paired-end sample: $SAMPLE_NAME"

    # Run Trim Galore! for paired-end data
    trim_galore --gzip \
        --path_to_cutadapt "$(which cutadapt)" \
        --cores 16 \
        --phred33 \
        --output_dir "$OUTPUT_DIR" \
        --paired \
        "$R1_FILE" "$R2_FILE"
done

