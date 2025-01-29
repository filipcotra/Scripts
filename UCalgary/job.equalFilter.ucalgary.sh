#!/bin/bash
#SBATCH --time=00-15:00
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.equalFilter.ucalgary
#SBATCH --output=./job.equalFilter.ucalgary.out
#SBATCH --error=./job.equalFilter.ucalgary.err
#SBATCH--mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

cd /home/cotra1/scratch/filteredMat1Variants 
 
module load python

strainArray=$(cat /home/cotra1/scratch/strainList.txt)

for strain in ${strainArray[@]}; do
         python equalFilter.py $strain.variants.vcf 
done
