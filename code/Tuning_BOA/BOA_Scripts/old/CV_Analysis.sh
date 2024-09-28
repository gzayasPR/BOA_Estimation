#!/bin/bash

# Set up the environment and variables
DIR=/blue/mateescu/gzayas97/Gabe_Thesis/2.Pop.Structure.GBC/Analysis
resultsBO=${DIR}/results/Breed.of.Origin/Tuning
my_bin=${DIR}/bin/Breed_of_Origin
NAME=UF.2024
F1=${DIR}/data/IDs/F1.ID

# Parameters to iterate over
folds=5
hidden_states=(2 5 10 15)
window_sizes=(5 10 15 20 30 40 50)
# Directory for storing the extracted ancestry data
comparison_dir=$resultsBO/CV_Results
mkdir -p "$comparison_dir"

# Loop through folds
for fold in $(seq -w 00 $((folds - 1))); do
    fold_dir="$comparison_dir/Fold_${fold}"
    mkdir -p "$fold_dir"
    awk '{print $1,$2}' "${resultsBO}/Mixed.Groups/Fold_${fold}/${NAME}.fam" > "$fold_dir/${NAME}.BO.ID"
    sed -i "1i PID IID" "$fold_dir/${NAME}.BO.ID"
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

        # Merge all runs for the chromosome in the fold side by side with headers
        merged_file="$fold_dir/merged_Chr${X}.txt"
        if [ ${#file_list[@]} -gt 0 ]; then
            paste "${file_list[@]}" > "$merged_file"
            headers=$(echo "${file_list[@]}" | sed -e 's/[^ ]*avg_ancestry_C[0-9][0-9]_//g'  -e 's/.txt//g' | tr ' ' '\t')
            sed -i "1i$headers" "$merged_file"
            
            # Combine the ID file with the merged ancestry file
            tmp_merged_file="${merged_file}_tmp"
            paste "$fold_dir/${NAME}.BO.ID" "$merged_file" > "$tmp_merged_file"
            mv "$tmp_merged_file" "$merged_file"

            # Filter the merged file for Angus and Brahman
            Angus_file="${fold_dir}/Angus_Chr${X}.txt"
            Brahman_file="${fold_dir}/Brahman_Chr${X}.txt"
            F1_file="${fold_dir}/F1_Chr${X}.txt"  # New line for F1 animals
            Angus_ids="${resultsBO}/Pure.Breed.Groups/Fold_${fold}/Angus_Fold_${fold}"
            Brahman_ids="${resultsBO}/Pure.Breed.Groups/Fold_${fold}/Brahman_Fold_${fold}"
            F1_ids="${F1}"  # New line for F1 IDs
            awk 'NR==FNR {id[$2]; next} $2 in id' "$Angus_ids" "$merged_file" > "$Angus_file"
            sed -i "1s/^/$(head -n1 "$merged_file")\n/" "$Angus_file"
            awk 'NR==FNR {id[$2]; next} $2 in id' "$Brahman_ids" "$merged_file" > "$Brahman_file"
            sed -i "1s/^/$(head -n1 "$merged_file")\n/" "$Brahman_file"
            awk 'NR==FNR {id[$2]; next} $2 in id' "$F1_ids" "$merged_file" > "$F1_file"  # New line for processing F1 data
            sed -i "1s/^/$(head -n1 "$merged_file")\n/" "$F1_file"  # New line to add header
        fi
    done
    rm ${fold_dir}/avg_ancestry*
done

echo "Ancestry data extracted and merged in $comparison_dir"


cd $comparison_dir


# File to store the parameter combinations
output_file="parameters.csv"

# Add header
echo "hidden_state,window_size" > "$output_file"

# Loop through hidden states and window sizes to generate combinations
for h_state in "${hidden_states[@]}"; do
    for w_size in "${window_sizes[@]}"; do
        echo "${h_state},${w_size}" >> "$output_file"
    done
done

echo "Parameter combinations saved to $output_file"


cd $comparison_dir
ml R
Rscript ${my_bin}/Tuning.scripts/CV_Analysis.r