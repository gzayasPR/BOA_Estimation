#!/bin/bash
#SBATCH --job-name=CHR_LAMPLD
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10gb
date; hostname; pwd

param_dir=$1
group=$2
X=$3
WINDOW_Size=$4
hidden_states=$5
proj_env=$6

source ${proj_env}
LAMPLD_Scripts=${LAMPLD_dir}
cd ${LAMPLD_Scripts}
mkdir -p ${param_dir}
cd ${param_dir}
cp ../chr.pos .
cp ../pop1.hap .
cp ../pop2.hap .
cp ../genofile.gen .
cd ${LAMPLD_Scripts}
touch ${param_dir}/LAMPLD${X}.out
perl run_LAMPLD2way.pl ${hidden_states} ${WINDOW_Size} ${param_dir}/chr.pos ${param_dir}/pop1.hap ${param_dir}/pop2.hap ${param_dir}/genofile.gen ${param_dir}/LAMPLD${X}.out
touch ${param_dir}/Chr.${X}.txt
perl convertLAMPLDout.pl ${param_dir}/LAMPLD${X}.out ${param_dir}/Chr.${X}.txt
touch ${param_dir}/LAMPLD.std_ancestry.txt
perl standardizeOutput.pl lamp-ld 2 ${param_dir}/Chr.${X}.txt ${param_dir}/LAMPLD.std_ancestry.txt
touch ${param_dir}/avg_ancestry.txt
perl averageAncestry.pl phased 2 ${param_dir}/LAMPLD.std_ancestry.txt ${param_dir}/avg_ancestry.txt
cd ${param_dir}
rm chr.pos 
rm pop1.hap 
rm pop2.hap
rm genofile.gen
echo "Chromosome ${X} finished"