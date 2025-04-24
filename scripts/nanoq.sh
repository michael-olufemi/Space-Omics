#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH -p cpu  # Partition
#SBATCH -t 20:00:00  # Job time limit
#SBATCH -o nanoq.log

source ~/.bashrc
conda activate nanoq


for file in *.fastq.gz; do
    nanoq -i "$file" -s -t 10 -vv > "${file%.fastq.gz}.nanoq"
done
