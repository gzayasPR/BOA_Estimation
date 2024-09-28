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
bash_out=$3
account_name=$4
folds=$5

# Display parameters for clarity
echo "Running pipeline with the following parameters:"
echo "      Project Environment: $proj_env"
echo "      Project Name: $NAME"
echo "      Bash Output Directory: $bash_out"
echo "      Account Name: $account_name"
echo "      Number of Folds: $folds"

# Source project environment and helper functions
source ${proj_env}
source ${my_code}/shell.functions.sh

# Set output directory for results
resultsBO=${my_results}/BOA_Tuning/

# Create directory for chromosome results
mkdir -p ${bash_out}/CHR

# Loop through folds
for fold in $(seq -w 00 $((folds - 1))); do
    # Initialize array to store job IDs
    job_ids=()
    fold_dir=${resultsBO}/Results/Fold_${fold}
    mkdir -p ${fold_dir}
    PurePath=${resultsBO}/Pure.Breed.Groups/Fold_$fold
    MixedPath=${resultsBO}/Mixed.Groups/Fold_$fold
    mkdir -p ${bash_out}/CHR/Fold_${fold}
    cd ${bash_out}/CHR/Fold_${fold}

    # Loop through chromosomes 1 to 29
    for X in $(seq -w 01 29); do
        chr_dir=${fold_dir}/${X}
        mkdir -p ${chr_dir}

        # Submit phasing job for the current chromosome
        job_id=$(sbatch --job-name="Phasing_${X}" \
            --output="Phasing_${X}.out" \
            --error="Phasing_${X}.err" \
            --account="${account_name}" \
            ${my_code}/Tuning_BOA/BOA_Scripts/CHR_Phasing.V3.sh ${chr_dir} ${MixedPath} ${PurePath} ${NAME} ${X} $proj_env | awk '{print $4}')

        # Store job ID
        job_ids+=($job_id)
    done

    # Wait for jobs in the current fold to complete
    echo "Waiting for jobs in fold $fold to complete..."
    for job_id in "${job_ids[@]}"; do
        wait_for_job_completion_less_output $job_id 5
    done

    echo "Fold $fold complete."
done

# Final wait to ensure all jobs are complete
echo "Waiting for all jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 10
done

echo "All jobs completed successfully!"
