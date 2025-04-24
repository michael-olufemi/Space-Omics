#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=1000G
#SBATCH -p cpu # Partition
#SBATCH --qos=long
#SBATCH -t 6-23:00:00  # Job time limit
#SBATCH -o download.log

python script.py
