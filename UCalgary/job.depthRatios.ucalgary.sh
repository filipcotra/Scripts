#!/bin/bash
#		d-hh:mm:ss
#SBATCH --time=1-00:00:00
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.depthRatios.ucalgary
#SBATCH --output=./depthRatios.out
#SBATCH --error=./depthRatios.err
#SBATCH--mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

cd .
 
module load samtools/1.12
module load python


rangeLower="1"
rangeUpper="1000"

strainList=$(cat /home/cotra1/scratch/RUFUSWorkflow/depth/SV_strainList.txt)
chromList=$(cat /home/cotra1/scratch/RUFUSWorkflow/depth/SV_chromList.txt)

chroLengthI="15072434"
chroLoopI="15072000"

chroLengthII="15279421"
chroLoopII="15279000"

chroLengthIII="13783801"
chroLoopIII="13783000"

chroLengthIV="17493829"
chroLoopIV="17493000"

chroLengthV="20924180"
chroLoopV="20924000"

chroLengthX="17718942"
chroLoopX="17718000"

for strain in ${strainList[@]}; do
	mkdir -p ${strain}Results
	echo "found strain"
	for chrom in ${chromList[@]}; do
		touch ${strain}Results/$strain.$chrom.depthRatio.txt
		loopBaseFormat=chroLoop$chrom
		lengthFormat=chroLength$chrom
		while [ $rangeUpper -le ${!loopBaseFormat} ]; do
			echo $rangeLower"-"$rangeUpper >> ${strain}Results/$strain.$chrom.depthRatio.txt
        		samtools depth -r $chrom:$rangeLower-$rangeUpper /home/cotra1/scratch/RUFUSWorkflow/bams/originalBAM/${strain}_trim_bwaMEM_sort_dedupped.bam |  awk 'BEGIN {sum = 0} {sum+=$3;} END { print "ave =",sum/NR;}'>> ${strain}Results/$strain.$chrom.depthRatio.txt
			rangeUpper=$(($rangeUpper+1000))
			rangeLower=$(($rangeLower+1000))
		done
		echo $rangeLower"-"${!lengthFormat} >> ${strain}Results/$strain.$chrom.depthRatio.txt
		samtools depth -r $chrom:$rangeLower-${!lengthFormat} /home/cotra1/scratch/RUFUSWorkflow/bams/originalBAM/${strain}_trim_bwaMEM_sort_dedupped.bam | awk 'BEGIN {sum = 0} {sum+=$3;} END { print "ave =",sum/NR;}' >> ${strain}Results/$strain.$chrom.depthRatio.txt
		rangeLower="1"
		rangeUpper="1000"

	#	python ratio.py ${strain}Results/$strain.$chrom.depthRatio.txt $(cat ${strain}Results/$strain.txt) ${strain}Results/$strain.$chrom.ratios.txt

	#	python normalization.py N2Results/N2.$chrom.ratios.txt ${strain}Results/$strain.$chrom.ratios.txt ${strain}Results/$strain.$chrom.normalizedDepth.txt

	done
done
