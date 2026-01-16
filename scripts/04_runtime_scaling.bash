#!/bin/bash
# ------------------------------------------------------------
# 04_runtime_scaling_array.bash
#   Runs split Stage 04 for each n in a SLURM array:
#     04a_simulate_geno.R
#     04b_compute_phi_amin.R
#     04c_run_rspectra.R
#     04d_run_rarpack.R
#     04e_combine_csv.R
# ------------------------------------------------------------

#SBATCH --job-name=runtime_scaling
#SBATCH --output=runtime_scaling_%A_%a.out
#SBATCH --array=1-10
#SBATCH --mem=250G
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:00
#SBATCH --mail-user=wx90@duke.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load R/4.4.3

# Values of n to benchmark (match array length)
N_VALUES=(10 20 50 100 200 500 1000 2000 5000 10000)

n=${N_VALUES[$SLURM_ARRAY_TASK_ID - 1]}

echo "SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID"
echo "Running n = $n"
echo "Working dir: $(pwd)"
echo "R: $(Rscript --version 2>&1)"

# Run split pipeline (each script is idempotent)
Rscript scripts/04a_simulate_geno.R $n
Rscript scripts/04b_compute_phi_amin.R $n
Rscript scripts/04c_run_rspectra.R $n
Rscript scripts/04d_run_rarpack.R $n
Rscript scripts/04e_combine_csv.R $n

echo "Finished n = $n"