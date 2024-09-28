#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20gb

# Log current date, hostname, and current working directory
date; hostname; pwd

# Input arguments
directory=$1    # Directory path
MixedPath=$2    # Path to mixed group files
PurePath=$3     # Path to purebred files
group=$4        # Group name
X=$5            # Chromosome number
proj_env=$6     # Project environment

# Source the project environment
source $proj_env

# Define path to R scripts
Rscripts=${my_code}/Breed_of_Origin/BOA_Scripts/Rscript/

# Create and navigate to the chromosome-specific directory
cd ${directory}
mkdir -p ./${X}
cd ${X}
echo "Processing group ${group}, chromosome ${X}"

# Load PLINK module for processing
module load plink/1.90b3

# Step 1: Run PLINK to recode the mixed group data for the specific chromosome
plink -cow -bfile ${MixedPath}/${group} -recode -chr ${X} --keep-allele-order --real-ref-alleles -out ${group}.${X}

# Step 2: Replace 'B' with 'T' in the .ped file
sed -i 's/\bB\b/T/g' ${group}.${X}.ped

# Step 3: Extract SNP IDs from the .map file
awk '{print $2}' ${group}.${X}.map > Chr${X}.snp

# Step 4: Recode purebred group data into VCF format for chromosome X
module load plink
plink -cow -bfile ${PurePath}/pop2 -chr ${X} --silent --keep-allele-order --real-ref-alleles -recode vcf --out pop2
plink -cow -bfile ${PurePath}/pop1 -chr ${X} --silent --keep-allele-order --real-ref-alleles -recode vcf --out pop1

# Step 5: Replace 'B' with 'T' in the VCF files
sed -i 's/\bB\b/T/g' pop1.vcf
sed -i 's/\bB\b/T/g' pop2.vcf

# Step 6: Run Beagle for phasing
module load beagle
beagle gt=pop2.vcf out=pop2.phased
beagle gt=pop1.vcf out=pop1.phased

# Unzip phased VCF files
gunzip pop2.phased.vcf.gz
gunzip pop1.phased.vcf.gz

# Step 7: Run R scripts to format haplotypes
module load R
Rscript ${Rscripts}/pop2.haps.format.R
Rscript ${Rscripts}/pop1.haps.format.R

# Clean up: Remove temporary phased VCF files
rm pop1.phased.vcf
rm pop2.phased.vcf

# Completion message
echo "Chromosome ${X} processing for group ${group} finished"
