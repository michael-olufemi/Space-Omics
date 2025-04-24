#!/bin/bash

# Base directory containing project directories
base_dir="/scratch3/workspace/michael_olufemi_student_uml_edu-HERRO/space_omics/RNA_analysis/infer_experiment"

# Threshold for determining strandedness
threshold=0.6

# Output summary file
summary_file="${base_dir}/strandedness_summary.txt"

# Initialize summary file
echo -e "Project_Directory\tSample_File\tStrandedness" > "$summary_file"

# Iterate over each project directory
for project_dir in "$base_dir"/*/; do
    project_name=$(basename "$project_dir")
    echo "Processing project: $project_name"

    # Iterate over each infer_expt.out file in the project directory
    for infer_file in "$project_dir"/*_infer_expt.out; do
        sample_name=$(basename "$infer_file" "_infer_expt.out")
        echo "  Analyzing sample: $sample_name"

        # Extract fractions from the infer_expt.out file
        frac_failed=$(awk '/Fraction of reads failed to determine:/ {print $NF}' "$infer_file")
        frac_1=$(awk '/Fraction of reads explained by "1\+\+,1--,2\+-\,2-\+":/ {print $NF}' "$infer_file")
        frac_2=$(awk '/Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":/ {print $NF}' "$infer_file")

        # Determine strandedness based on fractions
        if (( $(echo "$frac_1 > $threshold" | bc -l) )); then
            strandedness="reverse"
        elif (( $(echo "$frac_2 > $threshold" | bc -l) )); then
            strandedness="forward"
        else
            strandedness="none"
        fi

        # Append result to summary file
        echo -e "${project_name}\t${sample_name}\t${strandedness}" >> "$summary_file"
    done
done

echo "Strandedness analysis complete. Summary saved to $summary_file."
