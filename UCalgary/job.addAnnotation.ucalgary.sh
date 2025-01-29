#!/bin/bash
#SBATCH --time=00-15:00
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.annotationAdding.ucalgary
#SBATCH --output=./job.annotationAdding.ucalgary.out
#SBATCH --error=./job.annotationAdding.ucalgary.err
#SBATCH--mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

cd /home/cotra1/scratch/strainVariants

strainArray=$(cat /home/cotra1/scratch/strainList.txt)

for strain in ${strainArray[@]}; do
       python /home/cotra1/scratch/annotationAdder.py /home/cotra1/scratch/mat1_intragenicVariants.vcf $strain.variants.vcf /home/cotra1/scratch/annotatedStrainVariants/$strain.variants.vcf 
done
