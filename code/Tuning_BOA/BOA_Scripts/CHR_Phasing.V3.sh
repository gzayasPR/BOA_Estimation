#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
date;hostname;pwd
# 1 = directory
directory=$1
MixedPath=$2
PurePath=$3
group=$4
X=$5 
proj_env=$6

source $proj_env
Rscripts=${my_code}/Breed_of_Origin/BOA_Scripts/Rscript/

cd ${directory}
echo "${group} ${X}"
module load plink/1.90b3
plink -cow -bfile ${MixedPath}/${group}  -recode -chr ${X} --keep-allele-order --real-ref-alleles -out ${group}.${X}
sed -i 's/\bB\b/T/g' ${group}.${X}.ped
awk '{print $2}' ${group}.${X}.map > Chr${X}.snp
ml plink
plink -cow -bfile ${PurePath}/pop2 -chr ${X}  -silent  --keep-allele-order --real-ref-alleles -recode vcf  --out pop2
plink -cow -bfile ${PurePath}/pop1 -chr ${X}  -silent  --keep-allele-order --real-ref-alleles -recode vcf   --out pop1
sed -i 's/\bB\b/T/g' pop1.vcf
sed -i 's/\bB\b/T/g' pop2.vcf
ml beagle
beagle gt=pop2.vcf out=pop2.phased
beagle gt=pop1.vcf out=pop1.phased
gunzip pop2.phased.vcf.gz
gunzip pop1.phased.vcf.gz
ml R
Rscript ${Rscripts}/pop2.haps.format.R
Rscript ${Rscripts}/pop1.haps.format.R

rm pop1.phased.vcf
rm pop2.phased.vcf


cd ${LAMPLD_dir}
perl lait.pl lamp-ld 2 ${directory}/${group}.${X}.map ${directory}/${group}.${X}.ped ${directory}/Chr${X}.snp ${directory}/pop1.formatted_haplotypes.txt ${directory}/pop2.formatted_haplotypes.txt ${directory}