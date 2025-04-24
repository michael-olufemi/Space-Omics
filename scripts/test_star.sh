#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=500G
#SBATCH -p cpu  # Partition
#SBATCH --qos=long
#SBATCH -t 6-23:50:00  # Job time limit (4 days)


# Define directories
GENOME_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_genome"
OUTPUT_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_output/test_GLM-532"
LOG_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/logs"

# Create directories if they don't exist
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Define the input files
R1_FILE="./GLDS-532_rna-seq_GSM6594643_R1_raw_val_1.fq.gz"
R2_FILE="./GLDS-532_rna-seq_GSM6594643_R2_raw_val_2.fq.gz"

# Log file for debugging
LOG_FILE="${LOG_DIR}/star_test_GLM-532.log"

echo "Running STAR on test file: GLDS-532_rna-seq_GSM6594643" | tee -a "$LOG_FILE"
echo "Forward Read: $R1_FILE" | tee -a "$LOG_FILE"
echo "Reverse Read: $R2_FILE" | tee -a "$LOG_FILE"

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
     --outFileNamePrefix "${OUTPUT_DIR}/GLDS-532_" \
     --readFilesIn "$R1_FILE" "$R2_FILE" 2>&1 | tee -a "$LOG_FILE"

echo "STAR alignment completed for GLDS-532_rna-seq_GSM6594643" | tee -a "$LOG_FILE"
