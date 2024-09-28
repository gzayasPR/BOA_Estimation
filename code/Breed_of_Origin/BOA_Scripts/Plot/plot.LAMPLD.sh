#!/bin/bash
#SBATCH --account=mateescu 
#SBATCH --job-name=Plot_lamp
#SBATCH --output=Plot_lamp_%j.out
#SBATCH --error=errorfile_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=3999mb

cd /ufrc/mateescu/gzayas97/Thermo_Haplotypes/LAMPLD.PB/

mkdir Barplots
for X in {01..29};\
  do
  cd Barplots
  cp /ufrc/mateescu/gzayas97/Thermo_Haplotypes/LAMPLD.PB/$X/Chr.$X.txt .
  cp /ufrc/mateescu/gzayas97/Thermo_Haplotypes/LAMPLD.PB/$X/chr.pos chr$X.pos
  sed -e 's:\(.\):\1 :g' < Chr.$X.txt > $X.BO
  cd ../
  done
cd Barplots
module load R
Rscript /ufrc/mateescu/gzayas97/Thermo_Haplotypes/plotHAP.R
