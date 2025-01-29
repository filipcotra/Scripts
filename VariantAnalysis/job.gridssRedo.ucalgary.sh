#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-mtarailo
#SBATCH --job-name=./job.gridssRedo.ucalgary
#SBATCH --output=./job.gridssRedo.ucalgary.out
#SBATCH --error=./job.gridssRedo.ucalgary.err
#SBATCH --mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

cd /home/cotra1/scratch/gridss

module load bwa/0.7.17
module load r/4.1.0
module load samtools/1.12

/home/cotra1/projects/def-mtarailo/common/tools/gridss/gridss.sh --jar /home/cotra1/projects/def-mtarailo/common/tools/gridss/gridss-2.8.0-gridss-jar-with-dependencies.jar --reference /home/cotra1/projects/def-mtarailo/common/indexes/WS265_wormbase_july2020/c_elegens.PRJNA13758.WS265.genomic.fa --output CB4856.gridss.vcf.gz --assembly CB4856.gridss.assembly.bam --threads 4 --workingdir ./ --steps All --labels CB4856 /home/cotra1/projects/def-mtarailo/common/Hawaian.sorted.dedupped.bam
