#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=Run_BOA.Pipeline
#SBATCH --output=Run_BOA.Pipeline_%j.out
#SBATCH --error=Run_BOA.Pipeline_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb

# Log date, hostname, and working directory
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
NAME=UF.2024

# Define other parameters
states=10
WINDOW_Size=50
folds=3
missing_markers=0.01
missing_ind=0.01
afd=0.05


# Set up output directories
cd ${my_code}/Breed_of_Origin/
mkdir -p bash_out
bash_out=${my_code}/Breed_of_Origin/bash_out

# Step 1: Run PLINK Quality Control
mkdir -p ${bash_out}/Step1
cd ${bash_out}/Step1

job_ids=()
job_id=$(sbatch --job-name="${NAME}_Step1" \
                --output="${NAME}_Step1.out" \
                --error="${NAME}_Step1.err" \
                --account="${account_name}" \
                ${my_code}/Breed_of_Origin/BOA_Scripts/1.PLINK.QC.sh $proj_env $NAME $purelist $pop1_Pure $pop2_Pure $admix_list $genotypes $missing_markers $missing_ind $afd | awk '{print $4}')
job_ids+=($job_id)

# Wait for Step 1 to complete
echo "Waiting for Step 1 jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 5
done

# Check if required PLINK files were created
if [[ ! -f "${my_results}/Breed.of.Origin/Mixed.Groups/${NAME}.bed" || \
      ! -f "${my_results}/Breed.of.Origin/Pure.Breed.Groups/pop1.bed" || \
      ! -f "${my_results}/Breed.of.Origin/Pure.Breed.Groups/pop2.bed" ]]; then
    echo "Step 1 Error: Required PLINK files were not created."
    exit 1
fi

# Step 2: Run Phasing
mkdir -p ${bash_out}/Step2
cd ${bash_out}/Step2

job_id=$(sbatch --job-name="${NAME}_Step2" \
                --output="${NAME}_Step2.out" \
                --error="${NAME}_Step2.err" \
                --account="${account_name}" \
                ${my_code}/Breed_of_Origin/BOA_Scripts/2.run_Phasing_PB.sh $proj_env $NAME ${bash_out}/Step2 ${account_name} | awk '{print $4}')
job_ids+=($job_id)

# Wait for Step 2 to complete
echo "Waiting for Step 2 jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 5
done

# Check if required phasing files exist for each chromosome (01-29)
for X in $(seq -w 01 29); do
    if [[ ! -f "${my_results}/Breed.of.Origin/${NAME}/LAMPLD/${X}/${NAME}.${X}.map" || \
          ! -f "${my_results}/Breed.of.Origin/${NAME}/LAMPLD/${X}/${NAME}.${X}.ped" || \
          ! -f "${my_results}/Breed.of.Origin/${NAME}/LAMPLD/${X}/Chr${X}.snp" || \
          ! -f "${my_results}/Breed.of.Origin/${NAME}/LAMPLD/${X}/pop1.formatted_haplotypes.txt" || \
          ! -f "${my_results}/Breed.of.Origin/${NAME}/LAMPLD/${X}/pop2.formatted_haplotypes.txt" ]]; then
        echo "Error: Phasing step failed for chromosome ${X}. Required files were not created."
        exit 1
    fi
done

# Step 3: Run LAMPLD Analysis
mkdir -p ${bash_out}/Step3
cd ${bash_out}/Step3

job_id=$(sbatch --job-name="${NAME}_Step3" \
                --output="${NAME}_Step3.out" \
                --error="${NAME}_Step3.err" \
                --account="${account_name}" \
                ${my_code}/Breed_of_Origin/BOA_Scripts/3.run_LAMPLD.sh $proj_env $NAME ${states} ${WINDOW_Size} ${bash_out}/Step3 ${account_name} | awk '{print $4}')
job_ids+=($job_id)

# Wait for Step 3 to complete
echo "Waiting for Step 3 jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 5
done

# Step 4: Create PLINK files for analysis
mkdir -p ${bash_out}/Step4
cd ${bash_out}/Step4

job_id=$(sbatch --job-name="${NAME}_Step4" \
                --output="${NAME}_Step4.out" \
                --error="${NAME}_Step4.err" \
                --account="${account_name}" \
                ${my_code}/Breed_of_Origin/BOA_Scripts/4.Make.PLINK.sh $proj_env $NAME | awk '{print $4}')
job_ids+=($job_id)

# Wait for Step 4 to complete
echo "Waiting for Step 4 jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 5
done

# Completion message
echo "Pipeline completed successfully!"
