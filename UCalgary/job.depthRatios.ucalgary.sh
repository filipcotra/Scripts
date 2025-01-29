#!/bin/bash
#SBATCH --time=00-15:00
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.depthRatios.ucalgary
#SBATCH --output=./job.depthRatios.ucalgary.out
#SBATCH --error=./job.depthRatios.ucalgary.err
#SBATCH--mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

cd .
 
module load samtools/1.12

rangeLower="1"
rangeUpper="1000"

strainList=$(cat /home/cotra1/scratch/gridss/gridssList.txt)
chromList=$(cat /home/cotra1/scratch/chromList.txt)

chroLengthI="15072434"
chroLoopI="1507000"

chroLengthII="15279421"
chroLoopII="15279000"

chroLengthIII="13783801"
chroLoopIII="13783000"

chroLengthIV="17493829"
chroLoopIV="17493000"

chroLengthV="20924180"
chroLoopV="20924000"

for strain in ${strainList[@]}; do
	for chrom in ${chromList[@]}; do
		loopBaseFormat=chroLoop$chrom
		lengthFormat=chroLength$chrom
		while [ $rangeUpper -le ${!loopBaseFormat} ]; do
			echo $rangeLower"-"$rangeUpper >> $strain.$chrom.depthRatio.txt
        		samtools depth -r $chrom:$rangeLower-$rangeUpper /home/cotra1/scratch/gridss/$strain/$strain.bam |  awk 'BEGIN {sum = 0} {sum+=$3;} END { print "ave =",sum/NR;}'>> $strain.$chrom.depthRatio.txt
			rangeUpper=$(($rangeUpper+1000))
			rangeLower=$(($rangeLower+1000))
		done
		echo $rangeLower"-"${!lengthFormat} >> $strain.$chrom.depthRatio.txt
		samtools depth -r $chrom:$rangeLower-${!lengthFormat} /home/cotra1/scratch/gridss/$strain/$strain.bam | awk 'BEGIN {sum = 0} {sum+=$3;} END { print "ave =",sum/NR;}' >> $strain.$chrom.depthRatio.txt
		rangeLower="1"
		rangeUpper="1000"
	done
done
