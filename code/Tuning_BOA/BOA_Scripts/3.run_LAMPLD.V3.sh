#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb

# Log the current date, hostname, and working directory
date; hostname; pwd

# Input arguments
proj_env=$1
NAME=$2
hidden_states=($3)
window_sizes=($4)
folds=$5
bash_out=$6
account_name=$7

# Display parameters for clarity
echo "Running pipeline with the following parameters:"
echo "      Project Environment: $proj_env"
echo "      Project Name: $NAME"
echo "      Hidden States: ${hidden_states[*]}"
echo "      Window Sizes: ${window_sizes[*]}"
echo "      Number of Folds: $folds"
echo "      Bash Output Directory: $bash_out"
echo "      Account Name: $account_name"

# Source project environment and helper functions
source ${proj_env}
source ${my_code}/shell.functions.sh

# Set up result paths
resultsBO=${my_results}/BOA_Tuning/
PurePath=${resultsBO}/Pure.Breed.Groups
MixedPath=${resultsBO}/Mixed.Groups
job_ids=()

# Loop through folds
for fold in $(seq -w 00 $((folds - 1))); do
    fold_dir=${resultsBO}/Results/Fold_${fold}
    mkdir -p $bash_out/Fold_${fold}/
    cd  $bash_out/Fold_${fold}/
    
    # Loop through chromosomes 1 to 29
    for X in $(seq -w 01 29); do
        chr_dir=${fold_dir}/${X}
        mkdir -p $bash_out/Fold_${fold}/CHR
        cd  $bash_out/Fold_${fold}/CHR

        # Loop through hidden states and window sizes
        for h_state in "${hidden_states[@]}"; do
            for WINDOW_Size in "${window_sizes[@]}"; do
                param_dir=${chr_dir}/Window_${WINDOW_Size}_Hidden_${h_state}
                echo "Processing Fold: $fold, Chromosome: $X, Hidden States: $h_state, Window Size: $WINDOW_Size"

                # Submit job for current chromosome, window size, and hidden state
                job_id=$(sbatch --account="${account_name}" \
                    ${my_code}/Tuning_BOA/BOA_Scripts/CHR_LAMPLD.V3.sh  ${param_dir} $NAME $X $WINDOW_Size $h_state $proj_env | awk '{print $4}')
                
                # Store job ID
                job_ids+=($job_id)
            done
        done
    done

    # Wait for jobs in the current fold to complete
    for job_id in "${job_ids[@]}"; do
        wait_for_job_completion_less_output $job_id 5
    done
    echo "Fold ${fold} Complete"
done

# Final completion message
echo "All folds completed successfully!"
