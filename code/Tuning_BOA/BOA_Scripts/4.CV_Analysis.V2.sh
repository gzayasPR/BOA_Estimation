#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb

# Input arguments
proj_env=$1
NAME=$2
hidden_states=($3)
window_sizes=($4)
folds=$5
F1=$6

# Display parameters for clarity
echo "Running pipeline with the following parameters:"
echo "      Project Environment: $proj_env"
echo "      Project Name: $NAME"
echo "      Hidden States: ${hidden_states[*]}"
echo "      Window Sizes: ${window_sizes[*]}"
echo "      Number of Folds: $folds"
echo "      F1 File: $F1"

# Source project environment and helper functions
source ${proj_env}
source ${my_code}/shell.functions.sh

# Set up directories
resultsBO=${my_results}/BOA_Tuning/
comparison_dir=$resultsBO/CV_Results
mkdir -p "$comparison_dir"

# Loop through folds
for fold in $(seq -w 00 $((folds - 1))); do
    fold_dir="$comparison_dir/Fold_${fold}"
    mkdir -p "$fold_dir"

    # Generate BO.ID files
    awk '{print $1,$2}' "${resultsBO}/Mixed.Groups/Fold_${fold}/${NAME}.fam" > "$fold_dir/${NAME}.BO.ID"
    awk '{print $1,$2}' "${resultsBO}/Pure.Breed.Groups/Fold_${fold}/pop1_Fold_${fold}" > "$fold_dir/pop1.BO.ID"
    awk '{print $1,$2}' "${resultsBO}/Pure.Breed.Groups/Fold_${fold}/pop2_Fold_${fold}" > "$fold_dir/pop2.BO.ID"
    cp ${F1} "$fold_dir/F1.BO.ID"

    # Add headers to the files
    sed -i "1i PID IID" "$fold_dir/${NAME}.BO.ID"
    sed -i "1i PID IID" "$fold_dir/pop1.BO.ID"
    sed -i "1i PID IID" "$fold_dir/pop2.BO.ID"
    sed -i "1i PID IID" "$fold_dir/F1.BO.ID"

    # Loop through chromosomes
    for X in $(seq -w 1 29); do
        file_list=()

        # Loop through hidden states and window sizes
        for h_state in "${hidden_states[@]}"; do
            for WINDOW_Size in "${window_sizes[@]}"; do
                param_dir="${resultsBO}/Results/Fold_${fold}/${X}/Window_${WINDOW_Size}_Hidden_${h_state}"

                # Check if avg_ancestry.txt exists, then extract data
                ancestry_file="${param_dir}/avg_ancestry.txt"
                if [[ -f "$ancestry_file" ]]; then
                    output_file="$fold_dir/avg_ancestry_C${X}_W${WINDOW_Size}_H${h_state}.txt"
                    awk 'NR > 1 {print $2}' "$ancestry_file" > "$output_file"
                    file_list+=("$output_file")
                fi
            done
        done

        # Merge all runs for the chromosome in the fold
        merged_file="$fold_dir/merged_Chr${X}.txt"
        if [ ${#file_list[@]} -gt 0 ]; then
            paste "${file_list[@]}" > "$merged_file"
            headers=$(echo "${file_list[@]}" | sed -e 's/[^ ]*avg_ancestry_C[0-9][0-9]_//g' -e 's/.txt//g' | tr ' ' '\t')
            sed -i "1i$headers" "$merged_file"

            # Combine the ID file with the merged ancestry file
            tmp_merged_file="${merged_file}_tmp"
            paste "$fold_dir/${NAME}.BO.ID" "$merged_file" > "$tmp_merged_file"
            mv "$tmp_merged_file" "$merged_file"
        fi
    done

    # Clean up temporary files
    rm ${fold_dir}/avg_ancestry*
done

echo "Ancestry data extracted and merged in $comparison_dir"

# Change directory to comparison directory
cd $comparison_dir

# Create a CSV file to store parameter combinations
output_file="parameters.csv"
echo "hidden_state,window_size" > "$output_file"

# Loop through hidden states and window sizes to generate combinations
for h_state in "${hidden_states[@]}"; do
    for w_size in "${window_sizes[@]}"; do
        echo "${h_state},${w_size}" >> "$output_file"
    done
done

echo "Parameter combinations saved to $output_file"

# Run R script for analysis
ml R/4.3
Rscript ${my_code}/Tuning_BOA/BOA_Scripts/CV_Analysis.V2.R