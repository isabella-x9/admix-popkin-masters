# scripts/05_plot_scaling.R
# ------------------------------------------------------------
# Plot runtime scaling results from Stage 04
#
# Inputs:
#   output/runtime_scaling_n_<n>.csv
#   output/tmp/popkin_n_<n>_timing.csv
#
# Outputs:
#   output/runtime_scaling.pdf
#   output/popkin_timing_breakdown.pdf
#   output/runtime_scaling_with_popkin.pdf
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
COL_RSPECTRA <- "blue"
COL_RARPACK <- "red"
COL_RSPECTRA_OP <- "lightskyblue"
COL_EXPLICIT_TOTAL <- "pink2"

COL_POPKIN <- c(
  loadX  = "#4c78a8",
  popkin = "#f58518",
  popkinA = "#54a24b",
  write  = "#b279a2"
)

# ============================================================
# 1) Read runtime_scaling_n_<n>.csv
# ============================================================
rows <- lapply(N_VALUES, function(n) {
  f <- sprintf("output/runtime_scaling_n_%d.csv", n)
  d <- safe_read(f)
  
  if (is.null(d) || nrow(d) == 0) {
    return(data.frame(
      n = n,
      m = NA_real_,
      k = NA_real_,
      XTime = NA_real_,
      PhiTime = NA_real_,
      RSpectraExplicitEigs = NA_real_,
      RSpectraExplicitTotal = NA_real_,
      rARPACKOperatorEigs = NA_real_,
      rARPACKOperatorTotal = NA_real_,
      RSpectraOperatorEigs = NA_real_,
      RSpectraOperatorTotal = NA_real_,
      ExplicitMethodTotal = NA_real_,
      OperatorMethodTotal = NA_real_
    ))
  }
  
  data.frame(
    n = d$n[1],
    m = if ("m" %in% names(d)) d$m[1] else NA_real_,
    k = if ("k" %in% names(d)) d$k[1] else NA_real_,
    XTime = if ("XTime" %in% names(d)) d$XTime[1] else NA_real_,
    PhiTime = if ("PhiTime" %in% names(d)) d$PhiTime[1] else NA_real_,
    RSpectraExplicitEigs = if ("RSpectraExplicitEigs" %in% names(d)) d$RSpectraExplicitEigs[1] else NA_real_,
    RSpectraExplicitTotal = if ("RSpectraExplicitTotal" %in% names(d)) d$RSpectraExplicitTotal[1] else NA_real_,
    rARPACKOperatorEigs = if ("rARPACKOperatorEigs" %in% names(d)) d$rARPACKOperatorEigs[1] else NA_real_,
    rARPACKOperatorTotal = if ("rARPACKOperatorTotal" %in% names(d)) d$rARPACKOperatorTotal[1] else NA_real_,
    RSpectraOperatorEigs = if ("RSpectraOperatorEigs" %in% names(d)) d$RSpectraOperatorEigs[1] else NA_real_,
    RSpectraOperatorTotal = if ("RSpectraOperatorTotal" %in% names(d)) d$RSpectraOperatorTotal[1] else NA_real_,
    ExplicitMethodTotal = if ("ExplicitMethodTotal" %in% names(d)) d$ExplicitMethodTotal[1] else NA_real_,
    OperatorMethodTotal = if ("OperatorMethodTotal" %in% names(d)) d$OperatorMethodTotal[1] else NA_real_
  )
})

res <- do.call(rbind, rows)
res <- res[order(res$n), ]

# ============================================================
# 2) Runtime scaling plot: eigensolver totals only
# ============================================================
res_long <- data.frame(
  n = rep(res$n, 3),
  runtime_sec = c(
    res$RSpectraExplicitTotal,
    res$rARPACKOperatorTotal,
    res$RSpectraOperatorTotal
  ),
  method = factor(
    rep(
      c(
        "RSpectra (explicit Theta)",
        "rARPACK (operator Theta)",
        "RSpectra (operator Theta)"
      ),
      each = nrow(res)
    ),
    levels = c(
      "RSpectra (explicit Theta)",
      "rARPACK (operator Theta)",
      "RSpectra (operator Theta)"
    )
  )
)

res_long <- res_long[is.finite(res_long$runtime_sec) & res_long$runtime_sec > 0, , drop = FALSE]

if (nrow(res_long) > 0) {
  p_runtime <- ggplot(res_long, aes(x = n, y = runtime_sec, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c(
      "RSpectra (explicit Theta)" = COL_RSPECTRA,
      "rARPACK (operator Theta)" = COL_RARPACK,
      "RSpectra (operator Theta)" = COL_RSPECTRA_OP
    )) +
    labs(
      x = "Number of individuals (n)",
      y = "Runtime (seconds, log scale)",
      color = "Method"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave("output/runtime_scaling.pdf", p_runtime, width = 8, height = 5)
}

# ============================================================
# 3) Popkin timing breakdown plot
# ============================================================
pop_rows <- lapply(N_VALUES, function(n) {
  f <- sprintf("output/tmp/popkin_n_%d_timing.csv", n)
  d <- safe_read(f)
  
  if (is.null(d) || nrow(d) == 0) {
    return(data.frame(
      n = n,
      load_X_sec = NA_real_,
      popkin_sec = NA_real_,
      popkin_A_sec = NA_real_,
      write_outputs_sec = NA_real_,
      total_sec = NA_real_
    ))
  }
  
  data.frame(
    n = d$n[1],
    load_X_sec = if ("load_X_sec" %in% names(d)) d$load_X_sec[1] else NA_real_,
    popkin_sec = if ("popkin_sec" %in% names(d)) d$popkin_sec[1] else NA_real_,
    popkin_A_sec = if ("popkin_A_sec" %in% names(d)) d$popkin_A_sec[1] else NA_real_,
    write_outputs_sec = if ("write_outputs_sec" %in% names(d)) d$write_outputs_sec[1] else NA_real_,
    total_sec = if ("total_sec" %in% names(d)) d$total_sec[1] else NA_real_
  )
})

pop <- do.call(rbind, pop_rows)
pop <- pop[order(pop$n), ]

pop_long <- data.frame(
  n = rep(pop$n, 4),
  sec = c(
    pop$load_X_sec,
    pop$popkin_sec,
    pop$popkin_A_sec,
    pop$write_outputs_sec
  ),
  step = factor(
    rep(
      c(
        "Load X",
        "Popkin kinship",
        "Compute A_min",
        "Write outputs"
      ),
      each = nrow(pop)
    ),
    levels = c(
      "Load X",
      "Popkin kinship",
      "Compute A_min",
      "Write outputs"
    )
  )
)

pop_long <- pop_long[is.finite(pop_long$sec) & pop_long$sec > 0, , drop = FALSE]

if (nrow(pop_long) > 0) {
  p_pop <- ggplot(pop_long, aes(x = n, y = sec, color = step)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c(
      "Load X" = unname(COL_POPKIN["loadX"]),
      "Popkin kinship" = unname(COL_POPKIN["popkin"]),
      "Compute A_min" = unname(COL_POPKIN["popkinA"]),
      "Write outputs" = unname(COL_POPKIN["write"])
    )) +
    labs(
      title = "Runtime Breakdown of Popkin Kinship Estimation",
      x = "Number of individuals (n)",
      y = "Runtime (seconds, log scale)",
      color = "Step"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave("output/popkin_timing_breakdown.pdf", p_pop, width = 8, height = 5)
}

# ============================================================
# 4) Runtime scaling with Phi and explicit method total
# ============================================================
comp_long <- data.frame(
  n = rep(res$n, 5),
  runtime_sec = c(
    res$PhiTime,
    res$RSpectraExplicitTotal,
    res$rARPACKOperatorTotal,
    res$RSpectraOperatorTotal,
    res$ExplicitMethodTotal
  ),
  method = factor(
    rep(
      c(
        "Phi / A_min total",
        "RSpectra (explicit Theta)",
        "rARPACK (operator Theta)",
        "RSpectra (operator Theta)",
        "Explicit method total = Phi + RSpectra explicit"
      ),
      each = nrow(res)
    ),
    levels = c(
      "Phi / A_min total",
      "RSpectra (explicit Theta)",
      "rARPACK (operator Theta)",
      "RSpectra (operator Theta)",
      "Explicit method total = Phi + RSpectra explicit"
    )
  )
)

comp_long <- comp_long[is.finite(comp_long$runtime_sec) & comp_long$runtime_sec > 0, , drop = FALSE]

if (nrow(comp_long) > 0) {
  p_comp <- ggplot(comp_long, aes(x = n, y = runtime_sec, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c(
      "Phi / A_min total" = unname(COL_POPKIN["popkin"]),
      "RSpectra (explicit Theta)" = COL_RSPECTRA,
      "rARPACK (operator Theta)" = COL_RARPACK,
      "RSpectra (operator Theta)" = COL_RSPECTRA_OP,
      "Explicit method total = Phi + RSpectra explicit" = COL_EXPLICIT_TOTAL
    )) +
    labs(
      title = "Runtime Scaling with Phi",
      x = "Number of individuals (n)",
      y = "Runtime (seconds, log scale)",
      color = "Component / Method"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave("output/runtime_scaling_with_popkin.pdf", p_comp, width = 8.4, height = 5.4)
}

cat("Wrote: output/runtime_scaling.pdf\n")
cat("Wrote: output/popkin_timing_breakdown.pdf\n")
cat("Wrote: output/runtime_scaling_with_popkin.pdf\n")


