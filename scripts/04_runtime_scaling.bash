#!/bin/bash
# ------------------------------------------------------------
# 04_runtime_scaling_array.bash
#   - Run from login node:
#       bash scripts/04_runtime_scaling.bash submit
#   - SLURM runs this same file for each array task 
#   - Each step skips if its expected output exists
#   - Policy: compute explicit Phi + RSpectra up to n=50000,
#             but skip Phi/RSpectra at n=100000 (rARPACK only)
# ------------------------------------------------------------

set -euo pipefail
trap 'echo "[ERROR] Failed at line $LINENO. n=${n:-NA}. cmd=$BASH_COMMAND" >&2' ERR

module load R/4.4.3

# n values (array index is 1-based)
N_VALUES=(10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000)

# memory buckets by ARRAY INDEX (1-based)
SMALL_IDX="1-7"      # 10..1000
MED_IDX="8-10"       # 2000..10000
LARGE_IDX="11-12"    # 20000, 50000
XL_IDX="13"          # 100000

SMALL_MEM="32G"
MED_MEM="256G"
LARGE_MEM="256G"     # explicit Phi for 20k, 50k
XL_MEM="128G"        # operator-only for 100k

TIME_LIMIT="06:00:00"
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

  echo "[SUBMIT] XL array=${XL_IDX} mem=${XL_MEM}"
  sbatch \
    --job-name=rt_xl \
    --output="${LOG_DIR}/rt_xl_%A_%a.out" \
    --array="${XL_IDX}" \
    --mem="${XL_MEM}" \
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
  echo "Use: bash scripts/04_runtime_scaling_array.bash submit"
  exit 1
fi

idx=$((SLURM_ARRAY_TASK_ID - 1))
n="${N_VALUES[$idx]}"

# Skip explicit Phi/RSpectra ONLY at n=100000
PHI_CUTOFF=100000
DO_PHI=1
if [[ "${n}" -ge "${PHI_CUTOFF}" ]]; then
  DO_PHI=0
fi

echo "=============================================="
echo "[JOB ] ${SLURM_JOB_ID:-NA}  [TASK] ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] n=${n}  mem=${SLURM_MEM_PER_NODE:-NA}"
echo "[INFO] workdir=$(pwd)"
echo "[INFO] DO_PHI=${DO_PHI}"
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
if [[ "${DO_PHI}" -eq 0 ]]; then
  echo "[SKIP] Phi/A_min disabled at n=${n} (stop at 100000)"
else
  if [[ -f "${phi_out}" && -f "${amin_out}" ]]; then
    echo "[SKIP] ${phi_out} and ${amin_out}"
  else
    Rscript scripts/04b_compute_phi_amin.R "${n}"
  fi
fi

echo "[STEP] 04c RSpectra eigs"
if [[ "${DO_PHI}" -eq 0 ]]; then
  echo "[SKIP] RSpectra requires explicit Phi (skipped at n=${n})"
else
  if [[ -f "${rs_rt_out}" ]]; then
    echo "[SKIP] ${rs_rt_out}"
  else
    Rscript scripts/04c_run_rspectra.R "${n}"
  fi
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
