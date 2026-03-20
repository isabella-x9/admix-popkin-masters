# scripts/06b_plot_rmsd_diag.R
# ------------------------------------------------------------
# 06b: Eigenvector agreement diagnostics
#   - Raw RMSD by rank k
#
# Reads:
#   output/tmp/rspectra_operator_theta_n_<n>_vectors.tsv
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

COL_EXPLICIT <- "blue"
COL_OPERATOR <- "red"
COL_GRID <- "grey85"

op_path <- sprintf("output/tmp/rspectra_operator_theta_n_%d_vectors.tsv", n)
ex_path <- sprintf("output/tmp/rspectra_theta_n_%d_vectors.tsv", n)

if (!file.exists(op_path)) stop("Missing: ", op_path)
if (!file.exists(ex_path)) stop("Missing: ", ex_path)

op <- as.matrix(read.table(op_path, header = FALSE))
ex <- as.matrix(read.table(ex_path, header = FALSE))

k <- min(ncol(op), ncol(ex))

# Sign-align explicit eigenvectors to operator eigenvectors
for (j in seq_len(k)) {
  if (sum(op[, j] * ex[, j]) < 0) ex[, j] <- -ex[, j]
}

rmsd <- function(x, y) sqrt(mean((x - y)^2))
rmsd_by_k <- sapply(seq_len(k), function(j) rmsd(op[, j], ex[, j]))

dir.create("output", showWarnings = FALSE, recursive = TRUE)

out_pdf <- sprintf("output/eigenvector_rmsd_by_rank_n_%d.pdf", n)

df <- data.frame(k = seq_len(k), rmsd = rmsd_by_k)
df <- df[is.finite(df$rmsd) & df$rmsd > 0, , drop = FALSE]

if (nrow(df) == 0) {
  warning("No positive finite RMSD values computed; skipping plot.")
} else {
  ylim <- range(df$rmsd, finite = TRUE)
  
  p <- ggplot(df, aes(x = k, y = rmsd)) +
    geom_line(linewidth = 0.8, color = COL_OPERATOR) +
    geom_point(shape = 21, size = 2.4, stroke = 0.6,
               color = COL_OPERATOR, fill = COL_EXPLICIT) +
    scale_y_log10(limits = ylim) +
    scale_x_continuous(breaks = pretty(df$k)) +
    labs(
      x = "Rank k (1 = largest)",
      y = "RMSD(v_operator, v_explicit)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = COL_GRID),
      panel.grid.minor = element_line(color = COL_GRID)
    )
  
  ggsave(out_pdf, plot = p, width = 7.2, height = 5.2)
}

cat("Using:\n")
cat(" - operator: ", op_path, "\n", sep = "")
cat(" - explicit: ", ex_path, "\n", sep = "")
cat("Wrote: ", out_pdf, "\n", sep = "")



