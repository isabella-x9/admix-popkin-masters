# scripts/05b_plot_totals_comparison.R
# ------------------------------------------------------------
# Plot end-to-end runtime comparison using RSpectra only
#
# Totals:
#   Explicit method total = PhiTime + RSpectraExplicitTotal
#   Operator method total = XTime + RSpectraOperatorTotal
#
# Inputs:
#   output/runtime_scaling_n_<n>.csv
#
# Outputs:
#   output/runtime_totals_comparison_rspectra.csv
#   output/runtime_totals_comparison_rspectra.pdf
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

dir.create("output", showWarnings = FALSE, recursive = TRUE)

safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
}

N_VALUES <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000)

# ---- colors ----
COL_EXPLICIT <- "magenta4"
COL_OPERATOR <- "orchid1"

# ============================================================
# Read combined runtime files
# ============================================================
rows <- lapply(N_VALUES, function(n) {
  f <- sprintf("output/runtime_scaling_n_%d.csv", n)
  d <- safe_read(f)
  
  if (is.null(d) || nrow(d) == 0) {
    return(data.frame(
      n = n,
      XTime = NA_real_,
      PhiTime = NA_real_,
      RSpectraExplicitTotal = NA_real_,
      RSpectraOperatorTotal = NA_real_,
      ExplicitMethodTotal = NA_real_,
      OperatorMethodTotal = NA_real_
    ))
  }
  
  data.frame(
    n = d$n[1],
    XTime = if ("XTime" %in% names(d)) d$XTime[1] else NA_real_,
    PhiTime = if ("PhiTime" %in% names(d)) d$PhiTime[1] else NA_real_,
    RSpectraExplicitTotal = if ("RSpectraExplicitTotal" %in% names(d)) d$RSpectraExplicitTotal[1] else NA_real_,
    RSpectraOperatorTotal = if ("RSpectraOperatorTotal" %in% names(d)) d$RSpectraOperatorTotal[1] else NA_real_,
    ExplicitMethodTotal = if ("ExplicitMethodTotal" %in% names(d)) d$ExplicitMethodTotal[1] else NA_real_,
    OperatorMethodTotal = if ("OperatorMethodTotal" %in% names(d)) d$OperatorMethodTotal[1] else NA_real_
  )
})

res <- do.call(rbind, rows)
res <- res[order(res$n), ]

# Save comparison table
write.csv(
  res[, c(
    "n", "XTime", "PhiTime",
    "RSpectraExplicitTotal",
    "RSpectraOperatorTotal",
    "ExplicitMethodTotal",
    "OperatorMethodTotal"
  )],
  "output/runtime_totals_comparison_rspectra.csv",
  row.names = FALSE
)

# ============================================================
# Long format for plotting
# ============================================================
df_long <- rbind(
  data.frame(n = res$n, method = "Explicit method = Phi + RSpectra explicit", runtime_sec = res$ExplicitMethodTotal),
  data.frame(n = res$n, method = "Operator method = X + RSpectra operator", runtime_sec = res$OperatorMethodTotal)
)

df_long <- df_long[is.finite(df_long$runtime_sec) & df_long$runtime_sec > 0, , drop = FALSE]

df_long$method <- factor(
  df_long$method,
  levels = c(
    "Explicit method = Phi + RSpectra explicit",
    "Operator method = X + RSpectra operator"
  )
)

# ============================================================
# Plot
# ============================================================
if (nrow(df_long) == 0) {
  warning("No finite RSpectra runtime totals found; skipping totals comparison plot.")
} else {
  p <- ggplot(df_long, aes(x = n, y = runtime_sec, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c(
      "Explicit method = Phi + RSpectra explicit" = COL_EXPLICIT,
      "Operator method = X + RSpectra operator" = COL_OPERATOR
    )) +
    labs(
      title = "End-to-End Runtime Comparison (RSpectra Only)",
      subtitle = "Explicit vs operator-based RSpectra",
      x = "Number of individuals (n)",
      y = "Runtime (seconds, log scale)",
      color = "Method"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave("output/runtime_totals_comparison_rspectra.pdf", p, width = 8.2, height = 5.2)
}

cat("Wrote: output/runtime_totals_comparison_rspectra.csv\n")
cat("Wrote: output/runtime_totals_comparison_rspectra.pdf\n")
cat("Wrote: output/runtime_totals_comparison_rspectra.pdf\n")