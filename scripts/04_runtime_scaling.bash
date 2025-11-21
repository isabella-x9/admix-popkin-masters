#!/bin/bash
##SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=test
#SBATCH --output=test.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=wx90@duke.edu
#SBATCH --mail-type=END,FAIL 


module load R/4.4.3


time Rscript 04_runtime_scaling.R


module purge