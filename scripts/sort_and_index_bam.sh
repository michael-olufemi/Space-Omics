#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=100G  # Memory allocation
#SBATCH -p cpu  # Partition
#SBATCH --qos=long  # Quality of service
#SBATCH -t 3-00:00:00  # Job runtime (3 days)
#SBATCH -o sort_index_bam_${OSD_NAME}.log  # Output log file (per directory)

# Load necessary modules
 module load samtools/1.19.2 

# Define directories (passed from sbatch)
BAM_DIR="${BASE_BAM_DIR}/${OSD_NAME}"
SORTED_BAM_DIR="${SORTED_BAM_BASE_DIR}/${OSD_NAME}"

# Create output directory if not exists
mkdir -p "$SORTED_BAM_DIR"

echo "Processing directory: $BAM_DIR"

# Process each BAM file in the OSD directory
for BAM_FILE in "$BAM_DIR"/*Aligned.sortedByCoord.out.bam; do
    # Extract sample name
    SAMPLE_NAME=$(basename "$BAM_FILE" | sed 's/_Aligned.sortedByCoord.out.bam//')

    echo "Sorting BAM file for sample: $SAMPLE_NAME"

    # Sort BAM file
    SORTED_BAM="$SORTED_BAM_DIR/${SAMPLE_NAME}_Aligned.sortedByCoord_sorted.out.bam"
    samtools sort -m 3G --threads 32 -o "$SORTED_BAM" "$BAM_FILE"

    echo "Sorting completed for $SAMPLE_NAME"

    # Index the sorted BAM file
    echo "Indexing BAM file for sample: $SAMPLE_NAME"
    samtools index -@ 32 "$SORTED_BAM"

    echo "Indexing completed for $SAMPLE_NAME"
done

echo "Sorting and indexing completed for directory: $OSD_NAME"
