#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=500G
#SBATCH -p cpu  # Partition
#SBATCH --qos=long
#SBATCH -t 6-23:50:00  # Job time limit (4 days)
#SBATCH -o /scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/logs/star_alignment_%x.log
#SBATCH --job-name=star_align_%x



# Define directories
GENOME_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_genome"
BASE_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/trimmed"
OUTPUT_BASE_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_output"
LOG_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/logs"

# Retrieve dataset name from environment variable
DATASET_NAME=${DATASET_NAME}

# Define dataset and output directories
DATASET_DIR="${BASE_DIR}/${DATASET_NAME}"
OUTPUT_DIR="${OUTPUT_BASE_DIR}/${DATASET_NAME}"
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Define log file path for this project
LOG_FILE="${LOG_DIR}/star_alignment_${DATASET_NAME}.log"

# Find all R1 FASTQ files in the current dataset directory
FASTQ_FILES=("$DATASET_DIR"/*_R1_raw_val_1.fq.gz)
NUM_FILES=${#FASTQ_FILES[@]}

if [ $NUM_FILES -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $DATASET_DIR. Skipping." | tee -a "$LOG_FILE"
    exit 1
fi

# Process each paired-end sample
for R1_FILE in "${FASTQ_FILES[@]}"; do
    # Derive the corresponding R2 filename
    R2_FILE="${R1_FILE/_R1_raw_val_1.fq.gz/_R2_raw_val_2.fq.gz}"

    # Check if the R2 file exists
    if [ ! -f "$R2_FILE" ]; then
        echo "ERROR: Missing R2 file for $R1_FILE. Skipping." | tee -a "$LOG_FILE"
        continue
    fi

    # Extract sample name from the R1 filename
    SAMPLE_NAME=$(basename "$R1_FILE" | sed 's/_R1_raw_val_1.fq.gz//')

    echo "Processing sample: $SAMPLE_NAME" | tee -a "$LOG_FILE"
    echo "R1 File: $R1_FILE" | tee -a "$LOG_FILE"
    echo "R2 File: $R2_FILE" | tee -a "$LOG_FILE"

    # Run STAR alignment
    /scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR-2.7.11b/source/STAR --twopassMode Basic \
         --limitBAMsortRAM 65000000000 \
         --genomeDir "$GENOME_DIR" \
         --outSAMunmapped Within \
         --outFilterType BySJout \
         --outSAMattributes NH HI AS NM MD MC \
         --outFilterMultimapNmax 20 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverReadLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --sjdbScore 1 \
         --readFilesCommand zcat \
         --runThreadN 32 \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts \
         --outSAMheaderHD @HD VN:1.4 SO:coordinate \
         --outFileNamePrefix "${OUTPUT_DIR}/${SAMPLE_NAME}_" \
         --readFilesIn "$R1_FILE" "$R2_FILE" | tee -a "$LOG_FILE"

    echo "STAR alignment completed for $SAMPLE_NAME." | tee -a "$LOG_FILE"
done
