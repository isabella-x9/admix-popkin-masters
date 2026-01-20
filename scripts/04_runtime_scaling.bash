#!/bin/bash
# ------------------------------------------------------------
# 04_runtime_scaling_array.bash
#   -   Run from login node to launch jobs with memory set by n 
#   -   Executed automatically by SLURM for each array task
#   -   If a job is killed, rerunning continues from the last completed step 
# ------------------------------------------------------------

set -euo pipefail
trap 'echo "[ERROR] Failed at line $LINENO. n=${n:-NA}. cmd=$BASH_COMMAND" >&2' ERR

module load R/4.4.3

# n values (array index 1..10)
N_VALUES=(10 20 50 100 200 500 1000 2000 5000 10000)

# memory buckets by ARRAY INDEX (1-based)
SMALL_IDX="1-6"   # 10..500
MED_IDX="7-9"     # 1000..5000
LARGE_IDX="10"    # 10000

SMALL_MEM="32G"
MED_MEM="128G"
LARGE_MEM="250G"
TIME_LIMIT="04:00:00"

MAIL_USER="wx90@duke.edu"

LOG_DIR="logs"
mkdir -p "${LOG_DIR}" output/tmp output

# ---------- submit mode ----------
if [[ "${1:-}" == "submit" ]]; then
  echo "[SUBMIT] SMALL array=${SMALL_IDX} mem=${SMALL_MEM}"
  sbatch \
    --job-name=rt_small \
    --output="${LOG_DIR}/rt_small_%A_%a.out" \
    --array="${SMALL_IDX}" \
    --mem="${SMALL_MEM}" \
    --ntasks-per-node=1 \
    --time="${TIME_LIMIT}" \
    --mail-user="${MAIL_USER}" \
    --mail-type=END,FAIL \
    "$0"

  echo "[SUBMIT] MED array=${MED_IDX} mem=${MED_MEM}"
  sbatch \
    --job-name=rt_med \
    --output="${LOG_DIR}/rt_med_%A_%a.out" \
    --array="${MED_IDX}" \
    --mem="${MED_MEM}" \
    --ntasks-per-node=1 \
    --time="${TIME_LIMIT}" \
    --mail-user="${MAIL_USER}" \
    --mail-type=END,FAIL \
    "$0"

  echo "[SUBMIT] LARGE array=${LARGE_IDX} mem=${LARGE_MEM}"
  sbatch \
    --job-name=rt_large \
    --output="${LOG_DIR}/rt_large_%A_%a.out" \
    --array="${LARGE_IDX}" \
    --mem="${LARGE_MEM}" \
    --ntasks-per-node=1 \
    --time="${TIME_LIMIT}" \
    --mail-user="${MAIL_USER}" \
    --mail-type=END,FAIL \
    "$0"

  echo "[SUBMIT] Done. Use: squeue -u \$USER"
  exit 0
fi

# ---------- worker mode ----------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "[ERROR] Not running as a SLURM array task."
  echo "Use: bash scripts/04_runtime_scaling.bash submit"
  exit 1
fi

idx=$((SLURM_ARRAY_TASK_ID - 1))
n="${N_VALUES[$idx]}"

echo "=============================================="
echo "[JOB ] ${SLURM_JOB_ID:-NA}  [TASK] ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] n=${n}  mem=${SLURM_MEM_PER_NODE:-NA}"
echo "[INFO] workdir=$(pwd)"
echo "=============================================="

# expected outputs for resume/skip
geno_out="output/tmp/geno_n_${n}.tsv"
phi_out="output/tmp/Phi_n_${n}.tsv"
amin_out="output/tmp/Phi_n_${n}_Amin.rds"
rs_rt_out="output/tmp/rspectra_n_${n}_runtime.csv"
ra_rt_out="output/tmp/rarpack_n_${n}_runtime.csv"
final_out="output/runtime_scaling_n_${n}.csv"

echo "[STEP] 04a simulate geno"
if [[ -f "${geno_out}" ]]; then
  echo "[SKIP] ${geno_out}"
else
  Rscript scripts/04a_simulate_geno.R "${n}"
fi

echo "[STEP] 04b compute Phi/A_min"
if [[ -f "${phi_out}" && -f "${amin_out}" ]]; then
  echo "[SKIP] ${phi_out} and ${amin_out}"
else
  Rscript scripts/04b_compute_phi_amin.R "${n}"
fi

echo "[STEP] 04c RSpectra eigs"
if [[ -f "${rs_rt_out}" ]]; then
  echo "[SKIP] ${rs_rt_out}"
else
  Rscript scripts/04c_run_rspectra.R "${n}"
fi

echo "[STEP] 04d rARPACK eigs"
if [[ -f "${ra_rt_out}" ]]; then
  echo "[SKIP] ${ra_rt_out}"
else
  Rscript scripts/04d_run_rarpack.R "${n}"
fi

echo "[STEP] 04e combine csv"
if [[ -f "${final_out}" ]]; then
  echo "[SKIP] ${final_out}"
else
  Rscript scripts/04e_combine_csv.R "${n}"
fi

echo "[DONE] Finished n=${n}"