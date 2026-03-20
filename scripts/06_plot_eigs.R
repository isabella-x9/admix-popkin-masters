# scripts/06_plot_eigs.R
# ------------------------------------------------------------
# Parity plots for RSpectra only:
#   1) Eigenvalues parity: x = RSpectra operator Theta, y = RSpectra explicit Theta
#   2) Eigenvectors parity: k = 1..K, one plot per page, sign-aligned, with y = x
#   3) Corr-by-rank summary: k = 1..K
#
# Reads:
#   output/tmp/rspectra_operator_theta_n_<n>_values.tsv, _vectors.tsv
#   output/tmp/rspectra_theta_n_<n>_values.tsv, _vectors.tsv
#
# Writes:
#   output/eigenvalues_parity_n_<n>.pdf
#   output/eigenvectors_parity_n_<n>.pdf
#   output/eigenvector_corr_by_rank_n_<n>.pdf
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

N_VALUES <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

COL_EXPLICIT <- "blue"
COL_OPERATOR <- "red"
COL_DIAG <- "grey25"
PT_ALPHA <- 0.25

read_vec <- function(path) as.numeric(scan(path, quiet = TRUE))
read_mat <- function(path) as.matrix(read.table(path, header = FALSE))

have_files_for <- function(n) {
  v_op <- sprintf("output/tmp/rspectra_operator_theta_n_%d_values.tsv", n)
  u_op <- sprintf("output/tmp/rspectra_operator_theta_n_%d_vectors.tsv", n)
  v_ex <- sprintf("output/tmp/rspectra_theta_n_%d_values.tsv", n)
  u_ex <- sprintf("output/tmp/rspectra_theta_n_%d_vectors.tsv", n)
  file.exists(v_op) && file.exists(u_op) && file.exists(v_ex) && file.exists(u_ex)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  n <- as.integer(args[1])
  if (!have_files_for(n)) {
    stop("Need BOTH rspectra_operator_theta and rspectra_theta eig files for n = ", n)
  }
} else {
  n <- NA_integer_
  for (cand in rev(N_VALUES)) {
    if (have_files_for(cand)) {
      n <- cand
      break
    }
  }
  if (!is.finite(n)) stop("No n found with both RSpectra operator and explicit eig outputs.")
}

op_vals <- read_vec(sprintf("output/tmp/rspectra_operator_theta_n_%d_values.tsv", n))
ex_vals <- read_vec(sprintf("output/tmp/rspectra_theta_n_%d_values.tsv", n))

op_vecs <- read_mat(sprintf("output/tmp/rspectra_operator_theta_n_%d_vectors.tsv", n))
ex_vecs <- read_mat(sprintf("output/tmp/rspectra_theta_n_%d_vectors.tsv", n))

k <- min(length(op_vals), length(ex_vals), ncol(op_vecs), ncol(ex_vecs))
op_vals <- op_vals[seq_len(k)]
ex_vals <- ex_vals[seq_len(k)]
op_vecs <- op_vecs[, seq_len(k), drop = FALSE]
ex_vecs <- ex_vecs[, seq_len(k), drop = FALSE]

for (j in seq_len(k)) {
  s <- sum(op_vecs[, j] * ex_vecs[, j])
  if (is.finite(s) && s < 0) ex_vecs[, j] <- -ex_vecs[, j]
}

# (1) Eigenvalues parity
df_vals <- data.frame(
  k = seq_len(k),
  operator = op_vals,
  explicit = ex_vals
)

df_vals <- df_vals[
  is.finite(df_vals$operator) &
    is.finite(df_vals$explicit) &
    df_vals$operator > 0 &
    df_vals$explicit > 0,
  , drop = FALSE
]

if (nrow(df_vals) == 0) {
  warning("No positive finite eigenvalues found; skipping eigenvalues parity plot.")
} else {
  xlim <- range(df_vals$operator, finite = TRUE)
  ylim <- range(df_vals$explicit, finite = TRUE)
  
  label_ranks <- intersect(c(1, 2, 3, 5, 10), df_vals$k)
  df_lab <- df_vals[df_vals$k %in% label_ranks, , drop = FALSE]
  
  p_vals <- ggplot(df_vals, aes(x = operator, y = explicit)) +
    geom_point(
      shape = 21,
      size = 2.8,
      stroke = 0.7,
      color = COL_OPERATOR,
      fill = COL_EXPLICIT,
      alpha = 0.9
    ) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1.0, linetype = 2, color = COL_DIAG) +
    geom_text(
      data = df_lab,
      aes(label = paste0("k=", k)),
      vjust = -0.6,
      size = 3.4
    ) +
    scale_x_log10(limits = xlim) +
    scale_y_log10(limits = ylim) +
    coord_equal() +
    labs(
      x = "Eigenvalue (Operator Theta)",
      y = "Eigenvalue (Explicit Theta)"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    filename = sprintf("output/eigenvalues_parity_n_%d.pdf", n),
    plot = p_vals,
    width = 6.8,
    height = 5.8
  )
}

# (2) Eigenvectors parity
out_vec_pdf <- sprintf("output/eigenvectors_parity_n_%d.pdf", n)
cors <- rep(NA_real_, k)

pdf(out_vec_pdf, width = 7.2, height = 6.6)
for (j in seq_len(k)) {
  xj <- op_vecs[, j]
  yj <- ex_vecs[, j]
  
  r <- suppressWarnings(cor(xj, yj, use = "complete.obs"))
  cors[j] <- r
  
  dfj <- data.frame(x = xj, y = yj)
  dfj <- dfj[is.finite(dfj$x) & is.finite(dfj$y), , drop = FALSE]
  
  if (nrow(dfj) == 0) next
  
  rng <- range(c(dfj$x, dfj$y), finite = TRUE)
  
  p_j <- ggplot(dfj, aes(x = x, y = y)) +
    geom_point(
      shape = 21,
      size = 1.6,
      stroke = 0.2,
      color = COL_OPERATOR,
      fill = COL_EXPLICIT,
      alpha = PT_ALPHA
    ) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1.0, linetype = 2, color = COL_DIAG) +
    coord_equal(xlim = rng, ylim = rng) +
    labs(
      title = sprintf("k = %d", j),
      x = "Operator entry",
      y = "Explicit entry"
    ) +
    theme_bw(base_size = 12)
  
  print(p_j)
}
dev.off()

# (3) Corr-by-rank summary
df_cor <- data.frame(k = seq_len(k), corr = cors)
df_cor <- df_cor[is.finite(df_cor$corr), , drop = FALSE]

if (nrow(df_cor) == 0) {
  warning("No correlations computed; skipping corr-by-rank plot.")
} else {
  y_min <- min(0.0, min(df_cor$corr, na.rm = TRUE))
  
  p_cor <- ggplot(df_cor, aes(x = k, y = corr)) +
    geom_hline(yintercept = 1, linewidth = 1.0, linetype = 2, color = COL_DIAG) +
    geom_line(linewidth = 0.8, color = COL_OPERATOR) +
    geom_point(shape = 21, size = 2.4, stroke = 0.6, color = COL_OPERATOR, fill = COL_EXPLICIT) +
    scale_y_continuous(limits = c(y_min, 1.0)) +
    scale_x_continuous(breaks = pretty(df_cor$k)) +
    labs(
      x = "Rank k (1 = largest)",
      y = "Correlation"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    filename = sprintf("output/eigenvector_corr_by_rank_n_%d.pdf", n),
    plot = p_cor,
    width = 6.9,
    height = 4.9
  )
}

cat(
  "Using:\n",
  sprintf(" - operator: output/tmp/rspectra_operator_theta_n_%d_[values|vectors].tsv\n", n),
  sprintf(" - explicit: output/tmp/rspectra_theta_n_%d_[values|vectors].tsv\n", n),
  "Wrote:\n",
  sprintf(" - output/eigenvalues_parity_n_%d.pdf\n", n),
  sprintf(" - output/eigenvectors_parity_n_%d.pdf\n", n),
  sprintf(" - output/eigenvector_corr_by_rank_n_%d.pdf\n", n),
  sep = ""
)



