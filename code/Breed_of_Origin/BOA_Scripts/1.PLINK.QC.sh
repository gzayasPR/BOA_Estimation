#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20gb

# Log current date, hostname, and working directory
date; hostname; pwd

# Input arguments
proj_env=$1
NAME=$2
purelist=$3
pop1_Pure=$4
pop2_Pure=$5
admix_list=$6
genotypes=$7
missing_markers=$8
missing_ind=$9
afd=${10}

# Source project environment
source $proj_env

# Location of R scripts
Rscripts=${my_code}/Breed_of_Origin/BOA_Scripts/Rscript/

# Set up directories
cd ${my_results}
mkdir -p ./Breed.of.Origin/{Pure.Breed.Groups,Mixed.Groups}

# Check if directory setup succeeded
if [[ ! -d "${my_results}/Breed.of.Origin" ]]; then
    echo "Error: Failed to create directory ${my_results}/Breed.of.Origin"
    exit 1
fi

# Step 1: Quality Control (QC) on Entire Population
cd ${my_results}/Breed.of.Origin/Mixed.Groups
ml plink/1.90b3.39
plink -cow -bfile ${genotypes} -keep ${admix_list} -geno $missing_markers -recode -write-snplist --keep-allele-order --not-chr 0,30-33 --real-ref-alleles --silent -out ${NAME}.1

# Check if PLINK output was successfully created
if [[ ! -f "${NAME}.1.snplist" ]]; then
    echo "Error: PLINK failed to create ${NAME}.1.snplist"
    exit 1
fi

# Step 2: Create Purebred Groups
cd ${my_results}/Breed.of.Origin/Pure.Breed.Groups
plink -cow -bfile ${genotypes} -keep ${purelist} -recode -extract ${my_results}/Breed.of.Origin/Mixed.Groups/${NAME}.1.snplist --geno $missing_markers --not-chr 0,30-33 --keep-allele-order --real-ref-alleles --silent -out Pure

# Check if Pure.bed was successfully created
if [[ ! -f "Pure.ped" ]]; then
    echo "Error: PLINK failed to create Pure.bed"
    exit 1
fi

# Step 3: Calculate Allele Frequency Differences (AFD)
plink -cow -file Pure --within ${purelist} --freq --make-bed --keep-allele-order --real-ref-alleles --geno $missing_ind --silent -out Pure

# Check if frequency and bed files were created
if [[ ! -f "Pure.frq.strat" || ! -f "Pure.bed" ]]; then
    echo "Error: PLINK failed to create Pure.frq.strat or Pure.bed"
    exit 1
fi

# Run R script to filter SNPs with AFD > ${afd}
ml R
Rscript ${Rscripts}/Diff.snplist.R $afd

# Check if Diff.snplist was created
if [[ ! -f "${my_results}/Breed.of.Origin/Pure.Breed.Groups/DIFF.snplist" ]]; then
    echo "Error: Diff.snplist was not created by the Rscript"
    exit 1
fi

# Step 4: Create New Test Genotypes with Extracted SNPs
cd ${my_results}/Breed.of.Origin/Mixed.Groups
plink -cow -file ${NAME}.1 -keep ${admix_list} --geno $missing_ind -extract ${my_results}/Breed.of.Origin/Pure.Breed.Groups/DIFF.snplist -recode -make-bed --keep-allele-order --real-ref-alleles --silent -out ${NAME}

# Check if the test genotype file was created
if [[ ! -f "${NAME}.bed" ]]; then
    echo "Error: PLINK failed to create test genotype file ${NAME}.bed"
    exit 1
fi

# Step 5: Extract SNPs for Purebreds and Separate Angus and Brahman Purebreds
cd ${my_results}/Breed.of.Origin/Pure.Breed.Groups
plink -cow -bfile ${genotypes} -keep ${purelist} --extract ${my_results}/Breed.of.Origin/Pure.Breed.Groups/DIFF.snplist -make-bed --keep-allele-order --real-ref-alleles --silent -out Pure.FINAL

# Check if Pure.FINAL.bed was created
if [[ ! -f "Pure.FINAL.bed" ]]; then
    echo "Error: PLINK failed to create Pure.FINAL.bed"
    exit 1
fi

# Step 6: Extract SNPs for Angus (pop1) and Brahman (pop2)
plink -cow -bfile Pure.FINAL -keep ${pop1_Pure} -make-bed --keep-allele-order --real-ref-alleles --silent -out pop1
plink -cow -bfile Pure.FINAL -keep ${pop2_Pure} -make-bed --keep-allele-order --real-ref-alleles --silent -out pop2

# Check if pop1 and pop2 bed files were created
if [[ ! -f "pop1.bed" || ! -f "pop2.bed" ]]; then
    echo "Error: PLINK failed to create pop1.bed or pop2.bed"
    exit 1
fi

# Completion message
echo "Script completed successfully!"
