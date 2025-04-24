#!/bin/bash
#SBATCH -c 32                     # Number of Cores per Task
#SBATCH --mem=500G                # Memory allocation
#SBATCH -p cpu               # Partition
#SBATCH -t 23:50:00             # Job time limit (5 days)
#SBATCH -o bismark_mapping2.log    # Output log file

# Load necessary modules (adjust versions as needed)
module load bowtie2/2.4.2+py3.8.12
module load samtools/1.9+py2.7.18

# Define input and output directories
input_dir="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/new_trimmed"
output_dir="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/mapping"
genome_dir="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/ref_genome"  # Path to bisulfite genome

# Path to Bismark tools
bismark_path="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/Bismark-0.24.2/"

# Create output directory if it doesn't exist


# Change to input directory where FASTQ files are located
cd $input_dir

# Loop through each FASTQ file (assuming paired-end reads)
for file1 in *_R1_val_1.fq.gz; do
    base=$(basename $file1 _R1_val_1.fq.gz)
    file2="${base}_R2_val_2.fq.gz"

    # Run Bismark mapping with specified options
    ${bismark_path}bismark --bowtie2 --genome $genome_dir \
            -1 $file1 -2 $file2 \
            -o $output_dir 

    # Change to output directory for deduplication
    cd $output_dir

module load samtools/1.9+py2.7.18

    # Deduplicate the BAM file
    ${bismark_path}deduplicate_bismark --bam -p ${base}_bismark_bt2_pe.bam

    # Create cov.gz files using Bismark methylation extractor
    ${bismark_path}bismark_methylation_extractor ${base}_bismark_bt2_pe.deduplicated.bam \
                                  --bedGraph --gzip --multicore 6

    # Return to the input directory for the next file
    cd $input_dir
done

echo "Bismark mapping and methylation extraction complete!"
