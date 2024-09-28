#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200gb

# Log current date, hostname, and working directory
date; hostname; pwd

# Input arguments
proj_env=$1
NAME=$2

# Display parameters for clarity
echo "Running script with the following parameters:"
echo "Project environment: $proj_env"
echo "Project name: $NAME"

# Source the project environment
source ${proj_env}

# Define paths
resultsBO=${my_results}/Breed.of.Origin
Rscripts=${my_code}/Breed_of_Origin/BOA_Scripts/Rscript/
PLINK_MAKE_SCRIPTS=${my_code}/Breed_of_Origin/BOA_Scripts/Cnv_BreedOrigin.PLINK_files/

# Set up output directories
cd ${resultsBO}/${NAME}
mkdir -p ./Breed_of_Origin.files/MK.ped

# Process each chromosome (01 to 29)
cd ${resultsBO}/${NAME}/Breed_of_Origin.files/MK.ped
for X in {01..29}; do
    echo "Processing chromosome ${X}..."
    
    # Combine the file operations using a pipeline and avoid intermediate files
    cp ${resultsBO}/${NAME}/LAMPLD/$X/Chr.$X.txt $X.hap

    # Replace 0 with Angus (A) and 1 with Brahman (B), transpose and combine haplotypes in fewer steps
    sed -e 's/./& /g' -e 's/0/A/g' -e 's/1/B/g' < $X.hap | \
    awk -f ${PLINK_MAKE_SCRIPTS}/transpose.awk.txt | \
    sed 's/\(.\) \(.\)/\1\2/g' | \
    awk -f ${PLINK_MAKE_SCRIPTS}/transpose.awk.txt | \
    sed 's/ //g' | \
    sed -e 's/\(.\)/\1 /g' > $X.geno.txt

    echo "Finished processing chromosome ${X}"
done

# Copy map file
cp ${resultsBO}/Mixed.Groups/${NAME}.map ${NAME}.BO.map

# Step 2: Grab IIDs from the .ped file and save into a new file
awk '{print $1,$2,$3,$4,$5,$6}' ${resultsBO}/Mixed.Groups/${NAME}.ped > ${NAME}.BO.ID

# Run R script
module load R
Rscript ${PLINK_MAKE_SCRIPTS}/test9_F.R

# Final file handling
cd ${resultsBO}/${NAME}/Breed_of_Origin.files
cp ${resultsBO}/${NAME}/Breed_of_Origin.files/MK.ped/${NAME}.BO.ped .
cp ${resultsBO}/${NAME}/Breed_of_Origin.files/MK.ped/${NAME}.BO.map .

# Convert into a binary file using PLINK
module load plink/1.90b3
plink -cow -file ${NAME}.BO -make-bed --keep-allele-order --real-ref-alleles -out ${NAME}.BO

echo "Finished ${NAME}.BO"
