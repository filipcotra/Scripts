#!/bin/bash
#		d-hh:mm:ss
#SBATCH --time=0-15:00:00
#SBATCH --account=def-mtarailo
#SBATCH --job-name=editBamChrom.ucalgary
#SBATCH --output=./editBAM.out
#SBATCH --error=./editBAM.err
#SBATCH--mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

cd .

module load samtools

fileList=$(cat bamList.txt)

for file in ${fileList[@]}; do
	samtools view -h originalBAM/$file | sed 's@\tV\t@\t5\t@g' | sed 's@SN:V@SN:5@g' | sed 's@\tIV\t@\t4\t@g' | sed 's@SN:IV@SN:4@g' | sed 's@\tIII\t@\t3\t@g' | sed 's@SN:III@SN:3@g' | sed 's@\tII\t@\t2\t@g'| sed 's@SN:II@SN:2@g' | sed 's@\tI\t@\t1\t@g' | sed 's@SN:I@SN:1@g' > alteredBAM/$file.tmp.sam
	head alteredBAM/$file.tmp.sam
	samtools view -b alteredBAM/$file.tmp.sam > alteredBAM/altered_$file
done
