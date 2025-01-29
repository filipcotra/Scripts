#!/bin/bash
#SBATCH --time=02:00:00		#time (DD-HH:MM)
#SBATCH --account=def-mtarailo
#SBATCH --job-name=job.filtr.ucalgary
#SBATCH --output=./job.filtr.ucalgary.out
#SBATCH --error=./job.filtr.ucalgary.err
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
module load intel/2018.3
module load bedtools/2.27.1
hostname  >&2


sample=$(basename ./calling/varscan/N2_S1_snp_varscan.vcf | cut -f1 -d"_");
perl /home/cotra1/projects/def-mtarailo/common/tamaroi/pipeline/filtration_v6_aug2019.pl -database /home/cotra1/projects/def-mtarailo/common/database/database_celegans_april2019.vcf -vcf ./calling/varscan/N2_S1_snp_varscan.vcf -keepoutput ./filtration/N2_S1_snp_varscan_keep.tmp -output ./filtration/N2_S1_snp_varscan_excluded.tmp -threshold 2 -sample ${sample} -vcfannot ./annotation/coovar/${sample}*snp_varscan/categorized-gvs.gvf -ortholist /home/cotra1/projects/def-mtarailo/common/indexes/ortholist2/master_2.csv -biomart /home/cotra1/projects/def-mtarailo/common/indexes/Biomart/mart_export.txt -exon /home/cotra1/projects/def-mtarailo/common/indexes/Biomart/mart_export_exon.bed -call varscan -annot coovar;

head -n1 ./filtration/N2_S1_snp_varscan_keep.tmp > ./filtration/N2_S1_snp_varscan_keep.tsv
sed 1d ./filtration/N2_S1_snp_varscan_keep.tmp | sort | uniq >> ./filtration/N2_S1_snp_varscan_keep.tsv
rm ./filtration/N2_S1_snp_varscan_keep.tmp

head -n1 ./filtration/N2_S1_snp_varscan_excluded.tmp > ./filtration/N2_S1_snp_varscan_excluded.tsv
sed 1d ./filtration/N2_S1_snp_varscan_excluded.tmp | sort | uniq >> ./filtration/N2_S1_snp_varscan_excluded.tsv
rm ./filtration/N2_S1_snp_varscan_excluded.tmp


#sample=$(basename ./calling/varscan/N2_S1_snp_varscan.vcf | cut -f1 -d"_");

		