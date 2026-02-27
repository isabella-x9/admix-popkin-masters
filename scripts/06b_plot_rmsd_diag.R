# scripts/06b_plot_rmsd_diag.R
# ------------------------------------------------------------
# 06b: Eigenvector agreement diagnostics
#   - Raw RMSD by rank k
#
# Reads:
#   output/tmp/rarpack_n_<n>_vectors.tsv
#   output/tmp/rspectra_theta_n_<n>_vectors.tsv
#
# Writes:
#   output/eigenvector_rmsd_by_rank_n_<n>.pdf
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/06b_plot_rmsd_diag.R <n>")

n <- as.integer(args[1]) 

# ---------- COLORS (match previous scripts) ----------
COL_RARPACK  <- "#1f77b4"
COL_RSPECTRA <- "#ff7f0e"
COL_GRID     <- "grey85"

ra_path <- sprintf("output/tmp/rarpack_n_%d_vectors.tsv", n)
rs_path <- sprintf("output/tmp/rspectra_theta_n_%d_vectors.tsv", n)

if (!file.exists(ra_path)) stop("Missing: ", ra_path)
if (!file.exists(rs_path)) stop("Missing: ", rs_path)

ra_path <- sprintf("output/tmp/rarpack_n_%d_vectors.tsv", n)
rs_path <- sprintf("output/tmp/rspectra_theta_n_%d_vectors.tsv", n)

if (!file.exists(ra_path)) stop("Missing: ", ra_path)
if (!file.exists(rs_path)) stop("Missing: ", rs_path)

ra <- as.matrix(read.table(ra_path, header = FALSE))
rs <- as.matrix(read.table(rs_path, header = FALSE))

k <- min(ncol(ra), ncol(rs))

# Sign-align RSpectra eigenvectors to rARPACK
for (j in 1:k) {
  if (sum(ra[, j] * rs[, j]) < 0) rs[, j] <- -rs[, j]
}

rmsd <- function(x, y) sqrt(mean((x - y)^2))
rmsd_by_k <- sapply(1:k, function(j) rmsd(ra[, j], rs[, j]))

dir.create("output", showWarnings = FALSE, recursive = TRUE)

out_pdf <- sprintf("output/eigenvector_rmsd_by_rank_n_%d.pdf", n)

df <- data.frame(k = 1:k, rmsd = rmsd_by_k)
df <- df[is.finite(df$rmsd) & df$rmsd > 0, , drop = FALSE]

if (nrow(df) == 0) {
  warning("No positive finite RMSD values computed; skipping plot.")
} else {
  ylim <- range(df$rmsd, finite = TRUE)
  
  p <- ggplot(df, aes(x = k, y = rmsd)) +
    geom_line(linewidth = 0.8, color = COL_RARPACK) +
    geom_point(shape = 21, size = 2.4, stroke = 0.6, color = COL_RARPACK, fill = COL_RSPECTRA) +
    scale_y_log10(limits = ylim) +
    scale_x_continuous(breaks = pretty(df$k)) +
    labs(
      title = sprintf("Raw RMSD by eigenvector rank (n = %d)", n),
      x = "Rank k (1 = largest)",
      y = "RMSD(v_rARPACK, v_RSpectra)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = COL_GRID),
      panel.grid.minor = element_line(color = COL_GRID)
    )
  
  ggsave(out_pdf, plot = p, width = 7.2, height = 5.2)
}

cat("Wrote: ", out_pdf, "\n", sep = "")
