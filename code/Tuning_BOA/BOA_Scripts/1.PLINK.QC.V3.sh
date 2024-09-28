#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20gb

# Logging date, hostname, and working directory
date; hostname; pwd

# Input Arguments
proj_env=$1
NAME=$2
purelist=$3
pop1_Pure=$4
pop2_Pure=$5
MixList=$6
F1=$7
genotypes=$8
fold=$9
missing_markers=${10}
missing_ind=${11}
afd=${12}

# Display parameters for clarity
echo "Running pipeline with the following parameters:"
echo "      Project Environment: $proj_env"
echo "      Project Name: $NAME"
echo "      Pure List: $purelist"
echo "      Pop1 Pure: $pop1_Pure"
echo "      Pop2 Pure: $pop2_Pure"
echo "      Mixed List: $MixList"
echo "      F1 Genotype File: $F1"
echo "      Genotype Data: $genotypes"
echo "      Folds: $fold"
echo "      Missing Markers Threshold: $missing_markers"
echo "      Missing Individual Threshold: $missing_ind"
echo "      Allele Frequency Difference (AFD): $afd"

# Source project environment
source $proj_env
Rscripts=${my_code}/Breed_of_Origin/BOA_Scripts/Rscript/

# Create necessary directories
mkdir -p ${my_results}/BOA_Tuning/{Pure.Breed.Groups,Mixed.Groups,Results}

# Check if directories were created successfully
if [[ ! -d "${my_results}/BOA_Tuning" || ! -d "${my_results}/BOA_Tuning/Pure.Breed.Groups" || ! -d "${my_results}/BOA_Tuning/Mixed.Groups" ]]; then
    echo "Error: Failed to create necessary directories."
    exit 1
fi

# Step 1: QC on Entire Population
echo "Step 1: Running QC on entire population..."
cd ${my_results}/BOA_Tuning/Mixed.Groups
ml plink/1.90b3
plink -cow -bfile ${genotypes} -keep ${MixList} -geno $missing_markers -recode -write-snplist --keep-allele-order --real-ref-alleles --silent -out ${NAME}.1

# Check if the snplist was created
if [[ ! -f "${NAME}.1.snplist" ]]; then
    echo "Error: PLINK failed to create ${NAME}.1.snplist. Please check the log file: ${NAME}.1.log for more details."
    exit 1
fi

# Step 2: Create Purebred Groups
echo "Step 2: Creating purebred groups..."
cd ${my_results}/BOA_Tuning/Pure.Breed.Groups
plink -cow -bfile ${genotypes} -keep ${purelist} -recode -extract ${my_results}/BOA_Tuning/Mixed.Groups/${NAME}.1.snplist -geno $missing_markers --not-chr 0,30-33 --keep-allele-order --real-ref-alleles --silent -out Pure

# Check if Pure.bed was created
if [[ ! -f "Pure.ped" ]]; then
    echo "Error: PLINK failed to create Pure.bed. Please check the log file: Pure.log for more details."
    exit 1
fi

# Step 3: Calculate Allele Frequency Differences (AFD)
echo "Step 3: Calculating allele frequency differences (AFD)..."
plink -cow -file Pure --within ${purelist} --freq --make-bed --keep-allele-order --real-ref-alleles --mind $missing_ind --silent -out Pure
ml R
Rscript ${Rscripts}/Diff.snplist.R $afd

# Check if Diff.snplist was created
if [[ ! -f "${my_results}/BOA_Tuning/Pure.Breed.Groups/DIFF.snplist" ]]; then
    echo "Error: Diff.snplist was not created by Rscript."
    exit 1
fi

# Step 4: Create New Test Genotypes
echo "Step 4: Creating new test genotypes..."
cd ${my_results}/BOA_Tuning/Mixed.Groups
plink -cow -file ${NAME}.1 -keep ${MixList} --mind $missing_ind -extract ${my_results}/BOA_Tuning/Pure.Breed.Groups/DIFF.snplist -recode -make-bed --keep-allele-order --real-ref-alleles --silent -out ${NAME}

# Check if the test genotype file was created
if [[ ! -f "${NAME}.bed" ]]; then
    echo "Error: PLINK failed to create ${NAME}.bed. Please check the log file: ${NAME}.log for more details."
    exit 1
fi

# Step 5: Extract SNPs for Purebreds and Separate pop1 and pop2
echo "Step 5: Extracting SNPs for purebreds and separating pop1 and pop2..."
cd ${my_results}/BOA_Tuning/Pure.Breed.Groups
plink -cow -bfile ${genotypes} -keep ${purelist} --extract ${my_results}/BOA_Tuning/Pure.Breed.Groups/DIFF.snplist -make-bed --keep-allele-order --real-ref-alleles -recode --silent -out Pure.FINAL

# Check if Pure.FINAL.bed was created
if [[ ! -f "Pure.FINAL.bed" ]]; then
    echo "Error: PLINK failed to create Pure.FINAL.bed. Please check the log file: Pure.FINAL.log for more details."
    exit 1
fi

# Separate pop1 and pop2 purebreds
plink -cow -bfile Pure.FINAL -keep ${pop1_Pure} -make-bed --keep-allele-order --real-ref-alleles --silent -out pop1
plink -cow -bfile Pure.FINAL -keep ${pop2_Pure} -make-bed --keep-allele-order --real-ref-alleles --silent -out pop2

# Check if pop1.bed and pop2.bed were created
if [[ ! -f "pop1.bed" || ! -f "pop2.bed" ]]; then
    echo "Error: PLINK failed to create pop1.bed or pop2.bed. Please check the log files: pop1.log and pop2.log for more details."
    exit 1
fi

# Step 6: Shuffle and Split pop1 and pop2 into Folds
echo "Step 6: Shuffling and splitting pop1 and pop2 into folds..."
cd ${my_results}/BOA_Tuning/Pure.Breed.Groups

# Shuffle and split pop1 and pop2 into fold files with zero-padded numbers (00, 01, 02, etc.)
shuf pop1.fam > pop1_shuffled.txt
split -n l/${fold} -d  pop1_shuffled.txt pop1_Fold_

shuf pop2.fam > pop2_shuffled.txt
split -n l/${fold} -d  pop2_shuffled.txt pop2_Fold_

# Check if folds were created
if [[ ! -f "pop1_Fold_00" || ! -f "pop2_Fold_00" ]]; then
    echo "Error: Failed to create fold files for pop1 or pop2."
    exit 1
fi

# Step 7: Create Folds and Process Genotypes
for f in $(seq -w 00 $((fold - 1))); do
    echo "Processing fold $f..."
    
    # Create directory for the fold and move the fold files
    mkdir -p ${my_results}/BOA_Tuning/Pure.Breed.Groups/Fold_${f}
    cd ${my_results}/BOA_Tuning/Pure.Breed.Groups/Fold_${f}
    
    # Move fold files
    mv ../pop1_Fold_${f} .
    mv ../pop2_Fold_${f} .
    
    # Combine pop1 and pop2 fold files
    cat pop2_Fold_${f} pop1_Fold_${f} > Train.Pure.txt
    
    # Combine with F1 to create Train.ADMIX.txt
    cat Train.Pure.txt ${F1} > Train.ADMIX.txt
    
    # Check if Train.ADMIX.txt was created
    if [[ ! -f "Train.ADMIX.txt" ]]; then
        echo "Error: Failed to create Train.ADMIX.txt for fold ${f}."
        exit 1
    fi

    # Create PLINK bed files for pop1 and pop2 after removing fold samples
    plink -cow -bfile ../pop2 --remove pop2_Fold_${f} -make-bed --silent -out pop2
    plink -cow -bfile ../pop1 --remove pop1_Fold_${f} -make-bed --silent -out pop1

    # Check if bed files were created
    if [[ ! -f "pop2.bed" || ! -f "pop1.bed" ]]; then
        echo "Error: PLINK failed to create pop1.bed or pop2.bed for fold ${f}. Please check the log files: pop1.log and pop2.log for more details."
        exit 1
    fi
    
    # Step 8: Create Mixed Group Genotypes
    echo "Creating PLINK files for mixed group for fold ${f}..."
    mkdir -p ${my_results}/BOA_Tuning/Mixed.Groups/Fold_${f}
    cd ${my_results}/BOA_Tuning/Mixed.Groups/Fold_${f}
    
    plink -cow -bfile ../${NAME} --keep ${my_results}/BOA_Tuning/Pure.Breed.Groups/Fold_${f}/Train.ADMIX.txt -make-bed --silent -out ${NAME}

    # Check if PLINK file was created for the mixed group
    if [[ ! -f "${NAME}.bed" ]]; then
        echo "Error: PLINK failed to create ${NAME}.bed for fold ${f}. Please check the log file: ${NAME}.log for more details."
        exit 1
    fi
done

# Completion message
echo "Script completed successfully!"