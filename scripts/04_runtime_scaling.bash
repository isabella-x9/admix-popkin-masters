#!/bin/bash
# ------------------------------------------------------------
# 04_runtime_scaling.bash
#   - Run from login node:
#       bash scripts/04_runtime_scaling.bash submit
#   - SLURM runs this same file for each array task
#   - Each step skips if its expected output exists
#
# Pipeline:
#   04a   simulate genotype matrix X
#   04aa  time loading X
#   04b   compute explicit Phi and A_min
#   04c   RSpectra on explicit Theta
#   04d   rARPACK on operator Theta
#   04f   RSpectra on operator Theta
#   04e   combine runtime CSV
#
# Policy:
#   - Run up to n=50,000 only
#   - n=50,000 gets its own heavier resource request
# ------------------------------------------------------------

set -euo pipefail
trap 'echo "[ERROR] Failed at line $LINENO. n=${n:-NA}. cmd=$BASH_COMMAND" >&2' ERR

module load R/4.4.3

# n values (array index is 1-based)
N_VALUES=(10 20 50 100 200 500 1000 2000 5000 10000 20000 50000)

# memory buckets by ARRAY INDEX (1-based)
SMALL_IDX="1-7"      # 10..1000
MED_IDX="8-10"       # 2000..10000
LARGE_IDX="11"       # 20000
HUGE_IDX="12"        # 50000

SMALL_MEM="32G"
MED_MEM="256G"
LARGE_MEM="256G"
HUGE_MEM="320G"

SMALL_TIME="06:00:00"
MED_TIME="06:00:00"
LARGE_TIME="08:00:00"
HUGE_TIME="12:00:00"

MAIL_USER="wx90@duke.edu"

LOG_DIR="logs"
mkdir -p "${LOG_DIR}" output/tmp output

# ---------- submit mode ----------
if [[ "${1:-}" == "submit" ]]; then
  echo "[SUBMIT] SMALL array=${SMALL_IDX} mem=${SMALL_MEM} time=${SMALL_TIME}"
  sbatch \
    --job-name=rt_small \
    --output="${LOG_DIR}/rt_small_%A_%a.out" \
    --array="${SMALL_IDX}" \
    --mem="${SMALL_MEM}" \
    --ntasks-per-node=1 \
    --time="${SMALL_TIME}" \
    --mail-user="${MAIL_USER}" \
    --mail-type=END,FAIL \
    "$0"

  echo "[SUBMIT] MED array=${MED_IDX} mem=${MED_MEM} time=${MED_TIME}"
  sbatch \
    --job-name=rt_med \
    --output="${LOG_DIR}/rt_med_%A_%a.out" \
    --array="${MED_IDX}" \
    --mem="${MED_MEM}" \
    --ntasks-per-node=1 \
    --time="${MED_TIME}" \
    --mail-user="${MAIL_USER}" \
    --mail-type=END,FAIL \
    "$0"

  echo "[SUBMIT] LARGE array=${LARGE_IDX} mem=${LARGE_MEM} time=${LARGE_TIME}"
  sbatch \
    --job-name=rt_large \
    --output="${LOG_DIR}/rt_large_%A_%a.out" \
    --array="${LARGE_IDX}" \
    --mem="${LARGE_MEM}" \
    --ntasks-per-node=1 \
    --time="${LARGE_TIME}" \
    --mail-user="${MAIL_USER}" \
    --mail-type=END,FAIL \
    "$0"

  echo "[SUBMIT] HUGE array=${HUGE_IDX} mem=${HUGE_MEM} time=${HUGE_TIME}"
  sbatch \
    --job-name=rt_huge \
    --output="${LOG_DIR}/rt_huge_%A_%a.out" \
    --array="${HUGE_IDX}" \
    --mem="${HUGE_MEM}" \
    --ntasks-per-node=1 \
    --time="${HUGE_TIME}" \
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
geno_out="output/tmp/geno_n_${n}.rds"
xtime_out="output/tmp/X_n_${n}_timing.csv"
phi_out="output/tmp/Phi_n_${n}.rds"
amin_out="output/tmp/Phi_n_${n}_Amin.rds"

# Stage 04c
rs_rt_out="output/tmp/rspectra_theta_n_${n}_runtime.csv"

# Stage 04d
ra_rt_out="output/tmp/rarpack_n_${n}_runtime.csv"

# Stage 04f
rsop_rt_out="output/tmp/rspectra_operator_theta_n_${n}_runtime.csv"

# Stage 04e
final_out="output/runtime_scaling_n_${n}.csv"

echo "[STEP] 04a simulate geno"
if [[ -f "${geno_out}" ]]; then
  echo "[SKIP] ${geno_out}"
else
  Rscript scripts/04a_simulate_geno.R "${n}"
fi

echo "[STEP] 04aa time X load"
if [[ -f "${xtime_out}" ]]; then
  echo "[SKIP] ${xtime_out}"
else
  Rscript scripts/04aa_time_X.R "${n}"
fi

echo "[STEP] 04b compute Phi/A_min"
if [[ -f "${phi_out}" && -f "${amin_out}" ]]; then
  echo "[SKIP] ${phi_out} and ${amin_out}"
else
  Rscript scripts/04b_compute_phi_amin.R "${n}"
fi

echo "[STEP] 04c RSpectra eigs (explicit Theta)"
if [[ -f "${rs_rt_out}" ]]; then
  echo "[SKIP] ${rs_rt_out}"
else
  Rscript scripts/04c_run_rspectra.R "${n}"
fi

echo "[STEP] 04d rARPACK eigs (operator Theta)"
if [[ -f "${ra_rt_out}" ]]; then
  echo "[SKIP] ${ra_rt_out}"
else
  Rscript scripts/04d_run_rarpack.R "${n}"
fi

echo "[STEP] 04f RSpectra eigs (operator Theta)"
if [[ -f "${rsop_rt_out}" ]]; then
  echo "[SKIP] ${rsop_rt_out}"
else
  Rscript scripts/04f_run_rspectra_operator.R "${n}"
fi

echo "[STEP] 04e combine csv"
if [[ -f "${final_out}" ]]; then
  echo "[SKIP] ${final_out}"
else
  Rscript scripts/04e_combine_csv.R "${n}"
fi

echo "[DONE] Finished n=${n}"