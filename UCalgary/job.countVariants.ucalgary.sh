#!/bin/bash
#SBATCH --time=00-15:00
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.variantCounting.ucalgary
#SBATCH --output=./job.variantCounting.ucalgary.out
#SBATCH --error=./job.variantCounting.ucalgary.err
#SBATCH--mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

cd /home/cotra1/scratch

strainArray=$(cat /home/cotra1/scratch/strainList.txt)

for strain in ${strainArray[@]}; do
       python variantCounter.py /home/cotra1/scratch/filteredMat1Variants/$strain.variants.vcf /home/cotra1/scratch/variantCounts/mat1Counts.txt 
done
