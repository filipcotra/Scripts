#!/bin/bash
#SBATCH --time=02:00:00		#time (DD-HH:MM)
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.annotation.ucalgary
#SBATCH --output=./job.annotation.ucalgary.out
#SBATCH --error=./job.annotation.ucalgary.err
#SBATCH --mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=1		#number of mpi processes
#SBATCH --mem-per-cpu=10G	#memory; default unit is megabytes
#SBATCH --cpus-per-task=1

cd .

cd scratch

java -jar /home/cotra1/projects/def-mtarailo/common/tools/snpEff/snpEff.jar donwload Caenordhabditis_elegans

