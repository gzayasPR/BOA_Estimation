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
date;hostname;pwd

source ../project_env.sh
proj_env=${my_code}/project_env.sh

account_name="mateescu"
source ${my_code}/shell.functions.sh
##### Parameters/Inputs######
date;hostname;pwd
purelist=${my_data}/IDs/Pure.ID
pop1_Pure=${my_data}/IDs/Angus.ID
pop2_Pure=${my_data}/IDs/Brahman.ID
genotypes=${my_data}/Genotypes/UF_250K
admix_list=${my_data}/IDs/Mixed.Group.ID
F1=${my_data}/IDs/F1.ID
NAME=UF.2024
# Define the list of numbers
hidden_states=(2 3 5 10 15)
window_sizes=(5 10 15 20 30 40 50)
folds=3
# Quality Control
missing_markers=0.01
missing_ind=0.01
afd=0.05

cd ${my_code}/Tuning_BOA/
mkdir -p bash_out
bash_out=${my_code}/Tuning_BOA/bash_out

mkdir -p ${bash_out}/Step1
cd ${bash_out}/Step1

job_ids=()
job_id=$(sbatch --job-name="Step1_${NAME}" \
            --output="Step1_${NAME}1.out" \
            --error="Step1_${NAME}1.err" \
            --account="${account_name}" \
            ${my_code}/Tuning_BOA/BOA_Scripts/1.PLINK.QC.V3.sh $proj_env $NAME $purelist $pop1_Pure $pop2_Pure $admix_list $F1 $genotypes $folds $missing_markers $missing_ind $afd  | awk '{print $4}')
job_ids+=($job_id)
echo "Waiting for matrix jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 1
done

mkdir -p ${bash_out}/Step2
cd ${bash_out}/Step2
job_id=$(sbatch --job-name="Step2_${NAME}" \
            --output="Step2_${NAME}2.out" \
            --error="Step2_${NAME}2.err" \
            --account="${account_name}" \
            ${my_code}/Tuning_BOA/BOA_Scripts/2.run_Phasing_PB.V3.sh $proj_env $NAME ${bash_out}/Step2 ${account_name} $folds | awk '{print $4}')
job_ids+=($job_id)
echo "Waiting for matrix jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 10
done

mkdir -p ${bash_out}/Step3
cd ${bash_out}/Step3
job_id=$(sbatch --job-name="Step3_${NAME}" \
            --output="Step3_${NAME}3.out" \
            --error="Step3_${NAME}3.err" \
            --account="${account_name}" \
            ${my_code}/Tuning_BOA/BOA_Scripts/3.run_LAMPLD.V3.sh $proj_env $NAME "${hidden_states[*]}" "${window_sizes[*]}" $folds ${bash_out}/Step3 ${account_name} | awk '{print $4}')
job_ids+=($job_id)
echo "Waiting for matrix jobs to complete..."
for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 40
done

mkdir -p ${bash_out}/Step4
cd ${bash_out}/Step4

job_id=$(sbatch --job-name="Step4_${NAME}" \
            --output="Step4_${NAME}.out" \
            --error="Step4_${NAME}.err" \
            --account="${account_name}" \
            ${my_code}/Tuning_BOA/BOA_Scripts/4.CV_Analysis.V2.sh $proj_env $NAME "${hidden_states[*]}" "${window_sizes[*]}" $folds $F1 | awk '{print $4}')
job_ids+=($job_id)


for job_id in "${job_ids[@]}"; do
    wait_for_job_completion $job_id 2
done
