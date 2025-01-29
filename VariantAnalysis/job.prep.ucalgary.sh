#!/bin/bash
#SBATCH --time=02:00:00		#time (DD-HH:MM)
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.prep.ucalgary
#SBATCH --output=./job.prep.ucalgary.out
#SBATCH --error=./job.prep.ucalgary.err
#SBATCH --mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=1		#number of mpi processes
#SBATCH --mem-per-cpu=1024M	#memory; default unit is megabytes
#SBATCH --cpus-per-task=1

cd .

date

module avail

module load nixpkgs/16.09
module load perl/5.22.4 

#module load java/1.8.0_121 
module load java/1.7.0_80

module load samtools 

hostname  >&2

if [ ! -d ./calling ]; then mkdir ./calling; fi
if [ ! -d ./calling/varscan ]; then mkdir ./calling/varscan; fi
if [ ! -d ./annotation ]; then mkdir ./annotation; fi
if [ ! -d ./annotation/coovar ]; then mkdir ./annotation/coovar; fi
if [ ! -d ./filtration ]; then mkdir ./filtration; fi
if [ ! -e /home/cotra1/projects/def-mtarailo/common/indexes/WS265_wormbase/c_elegans.PRJNA13758.WS265.genomic.fa.fai ]; then samtools faidx /home/cotra1/projects/def-mtarailo/common/indexes/WS265_wormbase/c_elegans.PRJNA13758.WS265.genomic.fa; fi

echo " Prep gff3 for CooVar";
dir_gff3=$(dirname /home/cotra1/projects/def-mtarailo/common/indexes/WS265_wormbase/c_elegans.PRJNA13758.WS265.annotations.gff3);
base_gff3=$(basename /home/cotra1/projects/def-mtarailo/common/indexes/WS265_wormbase/c_elegans.PRJNA13758.WS265.annotations.gff3 .gff3);

grep -v '##' /home/cotra1/projects/def-mtarailo/common/indexes/WS265_wormbase/c_elegans.PRJNA13758.WS265.annotations.gff3 | grep 'ID=CDS' | grep "Parent" >> ${dir_gff3}/${base_gff3}_coovar.gff3;

	