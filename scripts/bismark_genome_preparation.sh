#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=500G 
#SBATCH -p cpu  # Partition
#SBATCH -t 20:00:00  # Job time limit
#SBATCH -o genome_preparation.log

module load bowtie/2.4.5 


./Bismark-0.24.2/bismark_genome_preparation ref_genome/
