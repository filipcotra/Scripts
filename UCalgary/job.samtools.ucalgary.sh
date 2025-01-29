#!/bin/bash
#SBATCH --time=00-15:00
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.samtools.ucalgary
#SBATCH --output=./job.samtools.ucalgary.out
#SBATCH --error=./job.samtools.ucalgary.err
#SBATCH--mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

cd .

module load samtools/1.12

samtools depth /home/cotra1/scratch/gridss/AB1/AB1.bam |  awk '{AB1sum+=$3;} END { print "AB1wholeAve =",AB1sum/NR;}'
samtools depth /home/cotra1/scratch/gridss/CB4856/CB4856.bam |  awk '{CBsum+=$3;} END { print "CB4856wholeAve =",CBsum/NR;}'
samtools depth /home/cotra1/scratch/gridss/GXW1/GXW1.bam |  awk '{GXsum+=$3;} END { print "GXW1wholeAve =",GXsum/NR;}'
samtools depth /home/cotra1/scratch/gridss/KR314/KR314.bam |  awk '{KRsum+=$3;} END { print "KR314wholeAve =",KRsum/NR;}'
samtools depth /home/cotra1/scratch/gridss/N2/N2.bam |  awk '{N2sum+=$3;} END { print "N2wholeAve =",N2sum/NR;}'
