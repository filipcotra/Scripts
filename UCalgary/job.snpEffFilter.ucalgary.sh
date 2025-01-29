#!/bin/bash
#SBATCH --time=02:00:00		#time (DD-HH:MM)
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.snpEffFilter.ucalgary
#SBATCH --output=./snpEffFilter.ucalgary.out
#SBATCH --error=./snpEffFilter.ucalgary.err
#SBATCH --mail-user=filip.cotra@ucalgary.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=1		#number of mpi processes
#SBATCH --mem-per-cpu=1024M	#memory; default unit is megabytes
#SBATCH --cpus-per-task=1

module load python

# Input the file names here
sampleFile="MTG441.90annotated.vcf"
orthoFile="info/Orthology.tsv"
goFile="info/Go_terms_1.tsv"
houseFile="info/InHouse.tsv"

# Code is run in the format: python PATH/snpEffFilter_November2022.py sampleFile orthoFile goFile houseFile [OPTIONS]
# For file paths, use the variables above. Change the path to the filter itself as needed. The options available are:
# -low and -mod, to include LOW and MODIFIER variants respectively. These can be added in whatever order you decide on.
python snpEffFilter_November2022.py $sampleFile $orthoFile $goFile $houseFile | sort $7 | uniq > uniqFiltered_${sampleFile}.snpeff

# Insert header
sed -i '1 i\CHROM\tPOS\tREF\tALT\tTYPE\tPRED\tGENE ID\tCHANGE\tLOCUS NAME\tGENE NAME\tHUMAN GENES\tPROTEIN NAME\tGO TERM\tGO DESC\tOMIM' uniqFiltered_${sampleFile}.snpeff
sed -i '2d' uniqFiltered_${sampleFile}.snpeff
