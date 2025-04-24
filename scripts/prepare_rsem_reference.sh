#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=200G
#SBATCH -p cpu # Partition
#SBATCH -t 23:00:00  # Job time limit
#SBATCH -o rsem_prepare_reference.log

module load conda/latest 
conda activate rsem_env


# Define paths
GENOME_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_genome"
FASTA_FILE="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF_FILE="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/Mus_musculus.GRCm39.107.gtf"


# Create output directory if not exists
mkdir -p $GENOME_DIR

# Run RSEM prepare reference
rsem-prepare-reference --gtf $GTF_FILE \
    $FASTA_FILE \
    $GENOME_DIR/RSEM_ref
