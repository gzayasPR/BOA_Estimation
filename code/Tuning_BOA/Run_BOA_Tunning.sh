#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=Run_BOA.Tunning
#SBATCH --output=Run_BOA.Tunning_%j.out
#SBATCH --error=Run_BOA.Tunning_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb

# Log current date, hostname, and working directory
date; hostname; pwd

# Load project environment and shell functions
source ../project_env.sh
proj_env=${my_code}/project_env.sh
source ${my_code}/shell.functions.sh

# Set parameters/inputs
account_name="mateescu"
purelist=${my_data}/IDs/Pure.ID
pop1_Pure=${my_data}/IDs/Angus.ID
pop2_Pure=${my_data}/IDs/Brahman.ID
genotypes=${my_data}/Genotypes/UF_250K
admix_list=${my_data}/IDs/Mixed.Group.ID
F1=${my_data}/IDs/F1.ID
NAME=UF.2024

# Define parameters for analysis
#hidden_states=(2 3 5 10 15)
#window_sizes=(5 10 15 20 30 40 50)
hidden_states=(2 3 )
window_sizes=(5 10 )
folds=3
missing_markers=0.01
missing_ind=0.01
afd=0.05

# Create necessary directories
cd ${my_code}/Tuning_BOA/
mkdir -p bash_out
bash_out=${my_code}/Tuning_BOA/bash_out

# Step 1: Quality Control (QC) 
mkdir -p ${bash_out}/Step1
cd ${bash_out}/Step1
echo "Starting Step 1 PLINK Quality Control"
job_ids=()
job_id=$(sbatch --job-name="Step1_${NAME}" \
                --output="Step1_${NAME}1.out" \
                --error="Step1_${NAME}1.err" \
                --account="${account_name}" \
                ${my_code}/Tuning_BOA/BOA_Scripts/1.PLINK.QC.V3.sh $proj_env $NAME $purelist $pop1_Pure $pop2_Pure $admix_list $F1 $genotypes $folds $missing_markers $missing_ind $afd | awk '{print $4}')
job_ids+=($job_id)

# Wait for QC jobs to complete
echo "Waiting for QC jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 5
done

# Check if required phasing files exist for each fold
for fold in $(seq -w 00 $((folds - 1))); do
    if [[ ! -f "${my_results}/BOA_Tuning/Mixed.Groups/Fold_${fold}/${NAME}.bed" || \
          ! -f "${my_results}/BOA_Tuning/Pure.Breed.Groups/Fold_${fold}/pop1.bed" || \
          ! -f "${my_results}/BOA_Tuning/Pure.Breed.Groups/Fold_${fold}/pop2.bed" ]]; then
        echo "Step 1 Error: PLINK for fold ${fold}. Required files were not created."
        exit 1
    fi
done

# Step 2: Phasing
mkdir -p ${bash_out}/Step2
cd ${bash_out}/Step2
echo "Starting Step 2 Phasing Purebreds using beagle"
job_ids=()
job_id=$(sbatch --job-name="Step2_${NAME}" \
                --output="Step2_${NAME}2.out" \
                --error="Step2_${NAME}2.err" \
                --account="${account_name}" \
                ${my_code}/Tuning_BOA/BOA_Scripts/2.run_Phasing_PB.V3.sh $proj_env $NAME ${bash_out}/Step2 ${account_name} $folds | awk '{print $4}')
job_ids+=($job_id)

# Wait for phasing jobs to complete
echo "Waiting for phasing jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 10
done

# Check if required phasing files exist for each fold
for fold in $(seq -w 00 $((folds - 1))); do
    for X in $(seq -w 01 29); do
        if [[ ! -f "${my_results}/BOA_Tuning/Results/Fold_${fold}/${X}/chr.pos" || \
            ! -f "${my_results}/BOA_Tuning/Results/Fold_${fold}/${X}/pop1.hap" || \
            ! -f "${my_results}/BOA_Tuning/Results/Fold_${fold}/${X}/pop2.hap" || \
            ! -f "${my_results}/BOA_Tuning/Results/Fold_${fold}/${X}/genofile.gen" ]]; then
            echo "Error: Phasing step failed for fold ${fold}. Required files were not created."
            exit 1
        fi
    done
done

# Step 3: LAMPLD Analysis
mkdir -p ${bash_out}/Step3
cd ${bash_out}/Step3
echo "Starting Step 3 Phasing Purebreds using Running LAMP-LD using different parameters"
job_ids=()
job_id=$(sbatch --job-name="Step3_${NAME}" \
                --output="Step3_${NAME}3.out" \
                --error="Step3_${NAME}3.err" \
                --account="${account_name}" \
                ${my_code}/Tuning_BOA/BOA_Scripts/3.run_LAMPLD.V3.sh $proj_env $NAME "${hidden_states[*]}" "${window_sizes[*]}" $folds ${bash_out}/Step3 ${account_name} | awk '{print $4}')
job_ids+=($job_id)

# Wait for LAMPLD jobs to complete
echo "Waiting for LAMPLD jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 15
done

# Step 4: Cross-Validation (CV) Analysis
mkdir -p ${bash_out}/Step4
cd ${bash_out}/Step4
job_ids=()
echo "Starting Step 4 Gathering Metrics for LAMP-LD"
job_id=$(sbatch --job-name="Step4_${NAME}" \
                --output="Step4_${NAME}.out" \
                --error="Step4_${NAME}.err" \
                --account="${account_name}" \
                ${my_code}/Tuning_BOA/BOA_Scripts/4.CV_Analysis.V2.sh $proj_env $NAME "${hidden_states[*]}" "${window_sizes[*]}" $folds $F1 | awk '{print $4}')
job_ids+=($job_id)

# Wait for CV jobs to complete
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 2
done

# Completion message
echo "Pipeline completed successfully!"
