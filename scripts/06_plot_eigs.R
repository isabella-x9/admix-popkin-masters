# scripts/06_plot_eigs.R
# ------------------------------------------------------------
# Parity plots:
#   1) Eigenvalues parity: x=rARPACK (operator Theta), y=RSpectra (explicit Theta), log–log, with y=x
#   2) Eigenvectors parity: k=1..K, one plot per page, sign-aligned, with y=x
#   3) Corr-by-rank summary: k=1..K
#
# Reads:
#   output/tmp/rarpack_n_<n>_values.tsv,  _vectors.tsv
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

# ---------- METHOD COLORS (match 05) ----------
COL_RARPACK  <- "#1f77b4"
COL_RSPECTRA <- "#ff7f0e"
COL_DIAG     <- "grey25"

PT_ALPHA <- 0.25

read_vec <- function(path) as.numeric(scan(path, quiet = TRUE))
read_mat <- function(path) as.matrix(read.table(path, header = FALSE))

method_prefix <- function(method) {
  method <- match.arg(method, c("rarpack", "rspectra"))
  if (method == "rarpack") "rarpack" else "rspectra_theta"
}

have_files_for <- function(n, method = c("rarpack", "rspectra")) {
  method <- match.arg(method)
  pref <- method_prefix(method)
  v <- sprintf("output/tmp/%s_n_%d_values.tsv", pref, n)
  u <- sprintf("output/tmp/%s_n_%d_vectors.tsv", pref, n)
  file.exists(v) && file.exists(u)
}

# Choose n: allow override "Rscript scripts/06_plot_eigs.R 20000"
N_VALUES <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  n <- as.integer(args[1])
  if (!have_files_for(n, "rarpack") || !have_files_for(n, "rspectra")) {
    stop("Need BOTH rarpack and rspectra_theta eig files for n = ", n)
  }
} else {
  n <- NA_integer_
  for (cand in rev(N_VALUES)) {
    if (have_files_for(cand, "rarpack") && have_files_for(cand, "rspectra")) {
      n <- cand
      break
    }
  }
  if (!is.finite(n)) stop("No n found with both rarpack and rspectra_theta eig outputs.")
}

# Load
ra_vals <- read_vec(sprintf("output/tmp/rarpack_n_%d_values.tsv", n))
rs_vals <- read_vec(sprintf("output/tmp/rspectra_theta_n_%d_values.tsv", n))

ra_vecs <- read_mat(sprintf("output/tmp/rarpack_n_%d_vectors.tsv", n))
rs_vecs <- read_mat(sprintf("output/tmp/rspectra_theta_n_%d_vectors.tsv", n))

# Ensure same k
k <- min(length(ra_vals), length(rs_vals), ncol(ra_vecs), ncol(rs_vecs))
ra_vals <- ra_vals[seq_len(k)]
rs_vals <- rs_vals[seq_len(k)]
ra_vecs <- ra_vecs[, seq_len(k), drop = FALSE]
rs_vecs <- rs_vecs[, seq_len(k), drop = FALSE]

# Sign alignment
for (j in seq_len(k)) {
  s <- sum(ra_vecs[, j] * rs_vecs[, j])
  if (is.finite(s) && s < 0) rs_vecs[, j] <- -rs_vecs[, j]
}

# ------------------------------------------------------------
# (1) Eigenvalues parity (log–log)
# ------------------------------------------------------------
df_vals <- data.frame(
  k = seq_len(k),
  rarpack = ra_vals,
  rspectra = rs_vals
)

df_vals <- df_vals[is.finite(df_vals$rarpack) & is.finite(df_vals$rspectra) &
                     df_vals$rarpack > 0 & df_vals$rspectra > 0, , drop = FALSE]

if (nrow(df_vals) == 0) {
  warning("No positive finite eigenvalues found; skipping eigenvalues parity plot.")
} else {
  xlim <- range(df_vals$rarpack, finite = TRUE)
  ylim <- range(df_vals$rspectra, finite = TRUE)
  
  label_ranks <- intersect(c(1, 2, 3, 5, 10), df_vals$k)
  df_lab <- df_vals[df_vals$k %in% label_ranks, , drop = FALSE]
  
  p_vals <- ggplot(df_vals, aes(x = rarpack, y = rspectra)) +
    geom_point(
      shape = 21,
      size = 2.8,
      stroke = 0.7,
      color = COL_RARPACK,
      fill = COL_RSPECTRA,
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
      title = sprintf("Eigenvalue parity (log–log, n = %d)", n),
      x = "Eigenvalue (rARPACK, operator Theta)",
      y = "Eigenvalue (RSpectra, explicit Theta)"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    filename = sprintf("output/eigenvalues_parity_n_%d.pdf", n),
    plot = p_vals,
    width = 6.8,
    height = 5.8
  )
}

# ------------------------------------------------------------
# (2) Eigenvectors parity: one plot per page
# ------------------------------------------------------------
out_vec_pdf <- sprintf("output/eigenvectors_parity_n_%d.pdf", n)
cors <- rep(NA_real_, k)

pdf(out_vec_pdf, width = 7.2, height = 6.6)
for (j in seq_len(k)) {
  xj <- ra_vecs[, j]
  yj <- rs_vecs[, j]
  
  r <- suppressWarnings(cor(xj, yj, use = "complete.obs"))
  cors[j] <- r
  
  rng <- range(c(xj, yj), finite = TRUE)
  xlim <- rng
  ylim <- rng
  
  dfj <- data.frame(x = xj, y = yj)
  dfj <- dfj[is.finite(dfj$x) & is.finite(dfj$y), , drop = FALSE]
  
  p_j <- ggplot(dfj, aes(x = x, y = y)) +
    geom_point(
      shape = 21,
      size = 1.6,
      stroke = 0.2,
      color = COL_RARPACK,
      fill = COL_RSPECTRA,
      alpha = PT_ALPHA
    ) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1.0, linetype = 2, color = COL_DIAG) +
    coord_equal(xlim = xlim, ylim = ylim) +
    annotate(
      "text",
      x = xlim[1] + 0.02 * diff(xlim),
      y = ylim[2] - 0.05 * diff(ylim),
      label = sprintf("corr = %.4f", r),
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    labs(
      title = sprintf("Eigenvector parity (n = %d): k = %d", n, j),
      x = "rARPACK entry (operator Theta)",
      y = "RSpectra entry (explicit Theta)"
    ) +
    theme_bw(base_size = 12)
  
  print(p_j)
}
dev.off()

# ------------------------------------------------------------
# (3) Corr-by-rank summary
# ------------------------------------------------------------
df_cor <- data.frame(k = seq_len(k), corr = cors)
df_cor <- df_cor[is.finite(df_cor$corr), , drop = FALSE]

if (nrow(df_cor) == 0) {
  warning("No correlations computed; skipping corr-by-rank plot.")
} else {
  y_min <- min(0.0, min(df_cor$corr, na.rm = TRUE))
  
  p_cor <- ggplot(df_cor, aes(x = k, y = corr)) +
    geom_hline(yintercept = 1, linewidth = 1.0, linetype = 2, color = COL_DIAG) +
    geom_line(linewidth = 0.8, color = COL_RARPACK) +
    geom_point(shape = 21, size = 2.4, stroke = 0.6, color = COL_RARPACK, fill = COL_RSPECTRA) +
    scale_y_continuous(limits = c(y_min, 1.0)) +
    scale_x_continuous(breaks = pretty(df_cor$k)) +
    labs(
      title = sprintf("Eigenvector agreement by rank (n = %d)", n),
      x = "Rank k (1 = largest)",
      y = "corr(rARPACK, RSpectra)"
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
  "Wrote:\n",
  sprintf(" - output/eigenvalues_parity_n_%d.pdf\n", n),
  sprintf(" - output/eigenvectors_parity_n_%d.pdf\n", n),
  sprintf(" - output/eigenvector_corr_by_rank_n_%d.pdf\n", n),
  sep = ""
)
