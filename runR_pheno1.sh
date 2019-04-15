#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=30000
#SBATCH --job-name=pheno1_GWAS
#SBATCH -o /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gmmat/pheno1/output/output.%j
#SBATCH -e /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gmmat/pheno1/output/errLog.%j
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1


cd $SLURM_SUBMIT_DIR
echo $PWD

module load r/3.5.2

echo "Rscript --vanilla GWAS.R pheno1/pheno1.gds pheno1/pheno1.txt /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gemma/gemma_output/pheno1.cXX.txt $1 pheno1/pheno1_results.txt"
Rscript --vanilla GWAS.R pheno1/pheno1.gds pheno1/pheno1.txt /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gemma/gemma_output/pheno1.cXX.txt $1 pheno1/pheno1_results.txt
