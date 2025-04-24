#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=600G
#SBATCH -p cpu # Partition
#SBATCH --qos=long
#SBATCH -t 6-23:00:00  # Job time limit
#SBATCH -o generate_star_ref.log

# Load the STAR module
module load star/2.7.11a

# Define variables
THREADS=16  # Number of threads to use
GENOME_DIR="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/STAR_genome"
FASTA_FILE="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF_FILE="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/Mus_musculus.GRCm39.107.gtf"
SJDB_OVERHANG=100  # For 101 bp reads, set to ReadLength - 1
RAM_LIMIT=550000000000  # 55 GB

# Create the genome directory if it doesn't exist
mkdir -p "$GENOME_DIR"

# Generate the STAR genome index
STAR --runThreadN "$THREADS" \
     --runMode genomeGenerate \
     --genomeDir "$GENOME_DIR" \
     --genomeFastaFiles "$FASTA_FILE" \
     --sjdbGTFfile "$GTF_FILE" \
     --sjdbOverhang "$SJDB_OVERHANG" \
     --limitGenomeGenerateRAM "$RAM_LIMIT"
