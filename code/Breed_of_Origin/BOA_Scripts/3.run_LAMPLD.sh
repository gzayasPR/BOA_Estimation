#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb

# Log current date, hostname, and working directory
date; hostname; pwd

# Input arguments
proj_env=$1
NAME=$2
states=$3
WINDOW_Size=$4
bash_out=$5
account_name=$6

# Display parameters for clarity
echo "Running script with the following parameters:"
echo "Project environment: $proj_env"
echo "Project name: $NAME"
echo "Number of hidden states: $states"
echo "Window size: $WINDOW_Size"
echo "Bash output directory: $bash_out"
echo "Account name: $account_name"

# Source the project environment and shell functions
source ${proj_env}
source ${my_code}/shell.functions.sh

# Define paths
resultsBO=${my_results}/Breed.of.Origin/
PurePath=${resultsBO}/Pure.Breed.Groups
MixedPath=${resultsBO}/Mixed.Groups
directory_input=${resultsBO}/${NAME}/LAMPLD

# Create necessary directories
mkdir -p ${resultsBO}/${NAME}
mkdir -p ${directory_input}
mkdir -p $bash_out/CHR

# Change to the bash output directory for chromosome-specific LAMPLD jobs
cd $bash_out/CHR

# Initialize an array to store job IDs
job_ids=()

# Submit jobs for chromosomes 01 to 29
for X in {01..29}; do
    echo "Processing chromosome ${X}"
    job_id=$(sbatch --job-name="LAMPLD_${X}" \
                    --output="LAMPLD_${X}.out" \
                    --error="LAMPLD_${X}.err" \
                    --account="${account_name}" \
                    ${my_code}/Breed_of_Origin/BOA_Scripts/CHR_LAMPLD.sh ${directory_input} ${MixedPath} ${PurePath} ${NAME} ${X} ${WINDOW_Size} ${states} $proj_env | awk '{print $4}')
    job_ids+=($job_id)
    sleep 2  # Small delay between job submissions
done

# Wait for all jobs to complete
echo "Waiting for matrix jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 5
done

echo "All jobs completed successfully!"
