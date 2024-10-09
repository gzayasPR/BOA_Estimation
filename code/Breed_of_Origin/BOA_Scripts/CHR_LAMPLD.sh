#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12gb

# Log the date, hostname, and current working directory
date; hostname; pwd

# Input arguments
directory=$1     # Directory path
MixedPath=$2     # Mixed group path
PurePath=$3      # Pure group path
group=$4         # Group name
X=$5             # Chromosome number
WINDOW_Size=$6   # Window size for LAMPLD
hidden_states=$7 # Number of hidden states
proj_env=$8      # Project environment

# Source the project environment
source ${proj_env}

# Define the path to LAMPLD scripts
LAMPLD_Scripts=${LAMPLD_dir}

# Change to the directory containing LAMPLD scripts
cd ${LAMPLD_Scripts}

# Step 1: Run LAMPLD mapping and SNP extraction
perl lait.pl lamp-ld 2 \
    ${directory}/${X}/${group}.${X}.map \
    ${directory}/${X}/${group}.${X}.ped \
    ${directory}/${X}/Chr${X}.snp \
    ${directory}/${X}/pop1.formatted_haplotypes.txt \
    ${directory}/${X}/pop2.formatted_haplotypes.txt \
    ${directory}/${X}

# Step 2: Run LAMPLD with 2-way ancestry analysis
perl run_LAMPLD2way.pl ${hidden_states} ${WINDOW_Size} \
    ${directory}/${X}/chr.pos \
    ${directory}/${X}/pop1.hap \
    ${directory}/${X}/pop2.hap \
    ${directory}/${X}/genofile.gen \
    ${directory}/${X}/LAMPLD${X}.out

# Step 3: Convert LAMPLD output to standard format
perl convertLAMPLDout.pl \
    ${directory}/${X}/LAMPLD${X}.out \
    ${directory}/${X}/Chr.${X}.txt

# Step 4: Standardize the output ancestry data
perl standardizeOutput.pl lamp-ld 2 \
    ${directory}/${X}/Chr.${X}.txt \
    ${directory}/${X}/LAMPLD.std_ancestry.txt

# Step 5: Calculate average ancestry across phased data
perl averageAncestry.pl phased 2 \
    ${directory}/${X}/LAMPLD.std_ancestry.txt \
    ${directory}/${X}/avg_ancestry.txt

# Print completion message for the chromosome
echo "Chromosome ${X} processing finished"