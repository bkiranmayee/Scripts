#!/bin/bash
 
#SBATCH -p assemble2
#SBATCH --mem 50000 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1


module load plink/2.00alM-2017-05-22

vcffile=$1
outfile=$2

plink2 --vcf $vcffile --const-fid --cow --export A-transpose --out $outfile

Rscript --vanilla vcfToAipl.R $outfile.traw $outfile

wait