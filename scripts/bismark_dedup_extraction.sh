#!/bin/bash
#SBATCH -c 32                     # Number of Cores per Task
#SBATCH --mem=500G                # Memory allocation
#SBATCH -p cpu-long                    # Partition
#SBATCH -t 6-20:00:00               # Job time limit (20 hours)
#SBATCH -o bismark_dedup_extraction.log   # Output log file

# Load necessary modules
module load samtools/1.9+py2.7.18

# Define input and output directories
input_dir="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/mapping"
output_dir="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/mapping"

# Path to Bismark tools
bismark_path="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/Bismark-0.24.2/"

# Change to the directory where BAM files are located
cd $input_dir

# Loop through each deduplicated BAM file and process
for bamfile in *_bismark_bt2_pe.bam; do
    base=$(basename $bamfile _bismark_bt2_pe.bam)

    # Deduplicate the BAM file
    ${bismark_path}deduplicate_bismark --bam ${bamfile}

    # Extract methylation information using Bismark methylation extractor
    ${bismark_path}bismark_methylation_extractor ${base}_bismark_bt2_pe.deduplicated.bam \
                                  --bedGraph --gzip --multicore 6

done

echo "Deduplication and methylation extraction complete!"

