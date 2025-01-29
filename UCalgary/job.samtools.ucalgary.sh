#!/bin/bash
#SBATCH --time=00-15:00
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.samtools.ucalgary
#SBATCH --output=./samtools.out
#SBATCH --error=./samtools.err
#SBATCH--mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

cd .

module load samtools/1.12

strainList=$(cat /home/cotra1/scratch/RUFUSWorkflow/depth/SV_strainList.txt)

for strain in ${strainList[@]}; do

	mkdir -p /home/cotra1/scratch/RUFUSWorkflow/depth/${strain}Results

	samtools depth /home/cotra1/scratch/RUFUSWorkflow/bams/originalBAM/${strain}_trim_bwaMEM_sort_dedupped.bam | awk '{sum+=$3;} END { print sum/NR;}' > /home/cotra1/scratch/RUFUSWorkflow/depth/${strain}Results/$strain.txt

done
