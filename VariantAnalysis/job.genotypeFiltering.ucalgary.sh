#!/bin/bash
#SBATCH --time=00-15:00
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.genotypeFiltering.ucalgary
#SBATCH --output=./job.genotypeFiltering.ucalgary.out
#SBATCH --error=./job.genotypeFiltering.ucalgary.err
#SBATCH--mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

cd /home/cotra1/scratch/cgh1VCFs 
 
module load vcftools

strainArray=$(cat /home/cotra1/scratch/strainList.txt)

for strain in ${strainArray[@]}; do
        vcftools --vcf /home/cotra1/scratch/cgh1_intragenicVariants.vcf --recode --recode-INFO-all --out $strain --indv $strain 
done
