#!/bin/bash

# Define input and output directories
input_dir="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/trimmed"
output_dir="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/mapping"
genome_dir="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/ref_genome"  # Path to bisulfite genome

# Path to Bismark tools
bismark_path="/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/Bismark-0.24.2/"


# Check if input directory is empty
if [ -z "$(ls -A $input_dir)" ]; then
    echo "Input directory is empty. No files to process."
    exit 1
fi

# Loop over all _R1_val_1.fq.gz files in the input directory
for file1 in ${input_dir}/*_R1_val_1.fq.gz; do
    # Get the base filename (removing the _R1_val_1.fq.gz suffix)
    base=$(basename $file1 _R1_val_1.fq.gz)
    
    # Define the corresponding _R2_val_2.fq.gz file
    file2="${input_dir}/${base}_R2_val_2.fq.gz"
    
    # Check if the _R2_val_2 file exists
    if [[ -f "$file2" ]]; then
        # Create a SLURM job script for this pair
        job_script="${output_dir}/${base}_bismark_mapping.sh"
        
        echo "#!/bin/bash" > $job_script
        echo "#SBATCH -c 32" >> $job_script                      # Number of cores per task
        echo "#SBATCH --mem=500G" >> $job_script                 # Memory allocation
        echo "#SBATCH -p cpu-long" >> $job_script                     # Partition
        echo "#SBATCH -t 6-23:00:00" >> $job_script                # Job time limit (20 hours)
        echo "#SBATCH -o ${output_dir}/${base}_bismark_mapping.log" >> $job_script  # Output log file

        # Load necessary modules
        echo "module load bowtie2/2.4.2+py3.8.12" >> $job_script
        echo "module load samtools/1.9+py2.7.18" >> $job_script

        # Run Bismark mapping with specified options
        echo "${bismark_path}bismark --bowtie2 --genome $genome_dir -1 $file1 -2 $file2 -o $output_dir" >> $job_script

        echo "echo 'Bismark mapping complete for $file1 and $file2'" >> $job_script
        
        # Submit the SLURM job
        job_id=$(sbatch $job_script | awk '{print $4}')
        
        if [[ $? -eq 0 ]]; then
            echo "Submitted batch job for $base with Job ID: $job_id"
        else
            echo "Error submitting batch job for $base"
        fi
    else
        echo "Warning: Corresponding R2 file for $file1 not found!"
    fi
done
