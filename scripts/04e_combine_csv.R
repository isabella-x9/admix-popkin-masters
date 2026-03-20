# scripts/04e_combine_csv.R
# ------------------------------------------------------------
# Combine per-n totals into output/runtime_scaling_n_<n>.csv
#
# Includes:
#   - X load time (04aa)
#   - popkin Phi build time (04b)
#   - RSpectra explicit Theta (04c)
#   - rARPACK operator Theta (04d)
#   - RSpectra operator Theta (04f)
#
# Derived totals:
#   - ExplicitMethodTotal = PhiTime + RSpectraExplicitTotal
#   - OperatorMethodTotal = XTime + RSpectraOperatorTotal
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04e_combine_csv.R <n>")

n <- as.integer(args[1])

geno_path <- sprintf("output/tmp/geno_n_%d.rds", n)
m <- NA_integer_
if (file.exists(geno_path)) {
  X <- readRDS(geno_path)
  m <- ncol(X)
}
k <- min(10, n - 1)

# Timing files
x_time <- sprintf("output/tmp/X_n_%d_timing.csv", n)    # 04aa
popkin_time <- sprintf("output/tmp/popkin_n_%d_timing.csv", n)    # 04b

rs_exp_time <- sprintf("output/tmp/rspectra_theta_n_%d_timing.csv", n)  # 04c
ra_op_time <- sprintf("output/tmp/rarpack_n_%d_timing.csv", n)    # 04d
rs_op_time <- sprintf("output/tmp/rspectra_operator_theta_n_%d_timing.csv", n) # 04f

# Fallback runtime files
rs_exp_rt <- sprintf("output/tmp/rspectra_theta_n_%d_runtime.csv", n)
ra_op_rt <- sprintf("output/tmp/rarpack_n_%d_runtime.csv", n)
rs_op_rt <- sprintf("output/tmp/rspectra_operator_theta_n_%d_runtime.csv", n)

# ---- Read times ----
XTime <- if (file.exists(x_time)) read.csv(x_time)$load_X_sec[1] else NA_real_

PhiTime <- if (file.exists(popkin_time)) {
  read.csv(popkin_time)$total_sec[1]
} else {
  NA_real_
}

RSpectraExplicitTotal <- if (file.exists(rs_exp_time)) {
  read.csv(rs_exp_time)$total_sec[1]
} else {
  if (file.exists(rs_exp_rt)) read.csv(rs_exp_rt)$runtime_sec[1] else NA_real_
}
RSpectraExplicitEigs <- if (file.exists(rs_exp_time)) {
  read.csv(rs_exp_time)$eigs_sec[1]
} else {
  NA_real_
}

rARPACKOperatorTotal <- if (file.exists(ra_op_time)) {
  read.csv(ra_op_time)$total_sec[1]
} else {
  if (file.exists(ra_op_rt)) read.csv(ra_op_rt)$runtime_sec[1] else NA_real_
}
rARPACKOperatorEigs <- if (file.exists(ra_op_time)) {
  read.csv(ra_op_time)$eigs_sec[1]
} else {
  NA_real_
}

RSpectraOperatorTotal <- if (file.exists(rs_op_time)) {
  read.csv(rs_op_time)$total_sec[1]
} else {
  if (file.exists(rs_op_rt)) read.csv(rs_op_rt)$runtime_sec[1] else NA_real_
}
RSpectraOperatorEigs <- if (file.exists(rs_op_time)) {
  read.csv(rs_op_time)$eigs_sec[1]
} else {
  NA_real_
}

# ---- Derived totals ----
ExplicitMethodTotal <- if (
  is.finite(PhiTime) &&
  is.finite(RSpectraExplicitTotal)
) {
  PhiTime + RSpectraExplicitTotal
} else {
  NA_real_
}

OperatorMethodTotal <- if (
  is.finite(XTime) &&
  is.finite(RSpectraOperatorTotal)
) {
  XTime + RSpectraOperatorTotal
} else {
  NA_real_
}

dir.create("output", showWarnings = FALSE, recursive = TRUE)
outfile <- sprintf("output/runtime_scaling_n_%d.csv", n)

write.csv(
  data.frame(
    n = n,
    m = m,
    k = k,
    
    XTime = XTime,
    PhiTime = PhiTime,
    
    RSpectraExplicitEigs = RSpectraExplicitEigs,
    RSpectraExplicitTotal = RSpectraExplicitTotal,
    
    rARPACKOperatorEigs = rARPACKOperatorEigs,
    rARPACKOperatorTotal = rARPACKOperatorTotal,
    
    RSpectraOperatorEigs = RSpectraOperatorEigs,
    RSpectraOperatorTotal = RSpectraOperatorTotal,
    
    ExplicitMethodTotal = ExplicitMethodTotal,
    OperatorMethodTotal = OperatorMethodTotal
  ),
  outfile,
  row.names = FALSE
)

cat("Wrote:", outfile, "\n")


