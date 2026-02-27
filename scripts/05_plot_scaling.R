# 05_plot_scaling.R
# ------------------------------------------------------------
# Plots:
#   1) Runtime scaling: RSpectra vs rARPACK
#   2) Popkin step breakdown (load X / popkin / popkin_A / write)
#   3) Load-time comparison: load X vs load Phi
#   4) NEW: Combined runtime plot (Popkin total + RSpectra total + rARPACK total)
#
# Reads:
#   output/runtime_scaling_n_<n>.csv  (created by 04e)
#   output/tmp/popkin_n_<n>_timing.csv
#   output/tmp/rspectra_n_<n>_timing.csv
#   output/tmp/rarpack_n_<n>_timing.csv
#
# Notes:
#   - All runtimes are in seconds
#   - base R log="xy" uses log10 axes (ticks at 10^k)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

N_VALUES <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ---------- GLOBAL COLORS ----------
COL_RARPACK  <- "#1f77b4"
COL_RSPECTRA <- "#ff7f0e"

COL_POPKIN <- c(
  loadX   = "#2ca02c",
  popkin  = "#d62728",
  popkinA = "#9467bd",
  write   = "#8c564b"
)

safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path)
}

# ---------- Runtime scaling ----------
rows <- lapply(N_VALUES, function(n) {
  f <- sprintf("output/runtime_scaling_n_%d.csv", n)
  d <- safe_read(f)
  if (is.null(d)) {
    return(data.frame(n = n, RSpectra = NA_real_, rARPACK = NA_real_))
  }
  data.frame(
    n = d$n[1],
    RSpectra = d$RSpectraTotal[1],
    rARPACK  = d$rARPACKTotal[1]
  )
})
res <- do.call(rbind, rows)

res_long <- rbind(
  data.frame(n = res$n, method = "rARPACK (operator)", runtime_sec = res$rARPACK),
  data.frame(n = res$n, method = "RSpectra (explicit Phi)", runtime_sec = res$RSpectra)
)

res_long <- res_long[is.finite(res_long$runtime_sec) & res_long$runtime_sec > 0 &
                       is.finite(res_long$n) & res_long$n > 0, , drop = FALSE]

if (nrow(res_long) == 0) {
  warning("No runtime scaling CSVs found in output/. Skipping runtime_scaling_full plot.")
} else {
  xlim_rt <- range(res_long$n, finite = TRUE)
  ylim_rt <- range(res_long$runtime_sec, finite = TRUE)
  
  p_runtime <- ggplot(res_long, aes(x = n, y = runtime_sec, color = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_log10(limits = xlim_rt) +
    scale_y_log10(limits = ylim_rt) +
    scale_color_manual(values = c(
      "rARPACK (operator)" = COL_RARPACK,
      "RSpectra (explicit Phi)" = COL_RSPECTRA
    )) +
    labs(
      title = "Runtime scaling (log10 axes): explicit Phi vs operator",
      x = "Sample size (n)",
      y = "Total runtime (seconds)",
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
  
  ggsave(
    filename = "output/runtime_scaling_full.pdf",
    plot = p_runtime,
    width = 7.2,
    height = 5.2
  )
}

# ---------- Popkin timing breakdown ----------
pop_rows <- lapply(N_VALUES, function(n) {
  f <- sprintf("output/tmp/popkin_n_%d_timing.csv", n)
  d <- safe_read(f)
  if (is.null(d)) return(NULL)
  d
})
pop <- do.call(rbind, pop_rows)

if (is.null(pop) || nrow(pop) == 0) {
  warning("No Popkin timing files found in output/tmp/. Skipping Popkin breakdown plot.")
} else {
  pop_long <- rbind(
    data.frame(n = pop$n, step = "load X", sec = pop$load_X_sec),
    data.frame(n = pop$n, step = "popkin()", sec = pop$popkin_sec),
    data.frame(n = pop$n, step = "popkin_A() (A_min)", sec = pop$popkin_A_sec),
    data.frame(n = pop$n, step = "write Phi / A_min", sec = pop$write_outputs_sec)
  )
  pop_long <- pop_long[is.finite(pop_long$sec) & pop_long$sec > 0 &
                         is.finite(pop_long$n) & pop_long$n > 0, , drop = FALSE]
  
  if (nrow(pop_long) == 0) {
    warning("Popkin timing data present but contains no positive finite values. Skipping Popkin breakdown plot.")
  } else {
    xlim_pop <- range(pop_long$n, finite = TRUE)
    ylim_pop <- range(pop_long$sec, finite = TRUE)
    
    p_pop <- ggplot(pop_long, aes(x = n, y = sec, color = step)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      scale_x_log10(limits = xlim_pop) +
      scale_y_log10(limits = ylim_pop) +
      #scale_color_manual(values = c(
      #  "load X" = COL_POPKIN["loadX"],
      #  "popkin()" = COL_POPKIN["popkin"],
      #  "popkin_A() (A_min)" = COL_POPKIN["popkinA"],
      #  "write Phi / A_min" = COL_POPKIN["write"]
      #)) +
      labs(
        title = "Popkin step breakdown (explicit Phi pipeline)",
        x = "Sample size (n)",
        y = "Seconds (log10 axes)",
        color = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "top")
    
    ggsave(
      filename = "output/popkin_timing_breakdown.pdf",
      plot = p_pop,
      width = 7.6,
      height = 5.2
    )
  }
}

# ---------- Load time: X vs Phi ----------
ra_rows <- lapply(N_VALUES, function(n) {
  d <- safe_read(sprintf("output/tmp/rarpack_n_%d_timing.csv", n))
  if (is.null(d)) return(NULL)
  data.frame(n = n, loadX = d$load_X_sec[1])
})
rs_rows <- lapply(N_VALUES, function(n) {
  d <- safe_read(sprintf("output/tmp/rspectra_n_%d_timing.csv", n))
  if (is.null(d)) return(NULL)
  data.frame(n = n, loadPhi = d$load_Phi_sec[1])
})

ra <- do.call(rbind, ra_rows)
rs <- do.call(rbind, rs_rows)

io_long <- rbind(
  data.frame(n = ra$n, what = "load X (rARPACK)", sec = ra$loadX),
  data.frame(n = rs$n, what = "load Phi (RSpectra)", sec = rs$loadPhi)
)

io_long <- io_long[is.finite(io_long$sec) & io_long$sec > 0 &
                     is.finite(io_long$n) & io_long$n > 0, , drop = FALSE]

if (nrow(io_long) == 0) {
  warning("No load-time timing files found. Skipping load_time_X_vs_Phi plot.")
} else {
  xlim_io <- range(io_long$n, finite = TRUE)
  ylim_io <- range(io_long$sec, finite = TRUE)
  
  p_io <- ggplot(io_long, aes(x = n, y = sec, color = what)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_log10(limits = xlim_io) +
    scale_y_log10(limits = ylim_io) +
    scale_color_manual(values = c(
      "load X (rARPACK)" = COL_RARPACK,
      "load Phi (RSpectra)" = COL_RSPECTRA
    )) +
    labs(
      title = "I/O cost (log10 axes): loading X vs loading Phi",
      x = "Sample size (n)",
      y = "Seconds (log10 axes)",
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
  
  ggsave(
    filename = "output/load_time_X_vs_Phi.pdf",
    plot = p_io,
    width = 7.6,
    height = 5.2
  )
}

# ---------- Combined runtime plot (Popkin total + RSpectra total + rARPACK total) ----------
# Popkin total runtime is pop$total_sec (explicit pipeline for building Phi/A_min)
# RSpectra/rARPACK totals are from output/runtime_scaling_n_<n>.csv (04e combine)
# All plotted in seconds, log="xy" in base R 

if (!is.null(pop) && nrow(pop) > 0) {
  pop_total <- data.frame(n = pop$n, PopkinTotal = pop$total_sec)
} else {
  pop_total <- data.frame(n = numeric(0), PopkinTotal = numeric(0))
}

combo <- merge(
  res[, c("n", "RSpectra", "rARPACK")],
  pop_total,
  by = "n",
  all = TRUE
)
combo <- combo[order(combo$n), ]

combo$ExplicitTotal <- NA_real_
okE <- is.finite(combo$PopkinTotal) & combo$PopkinTotal > 0 &
  is.finite(combo$RSpectra) & combo$RSpectra > 0
combo$ExplicitTotal[okE] <- combo$PopkinTotal[okE] + combo$RSpectra[okE]

ok_any <- (is.finite(combo$PopkinTotal)   & combo$PopkinTotal   > 0) |
  (is.finite(combo$RSpectra)      & combo$RSpectra      > 0) |
  (is.finite(combo$rARPACK)       & combo$rARPACK       > 0) |
  (is.finite(combo$ExplicitTotal) & combo$ExplicitTotal > 0)
combo <- combo[ok_any, ]

combo_long <- rbind(
  data.frame(n = combo$n, series = "Popkin total (build Phi / A_min)", sec = combo$PopkinTotal),
  data.frame(n = combo$n, series = "RSpectra (explicit Phi) total", sec = combo$RSpectra),
  data.frame(n = combo$n, series = "rARPACK (operator) total", sec = combo$rARPACK),
  data.frame(n = combo$n, series = "Explicit pipeline total = Popkin + RSpectra", sec = combo$ExplicitTotal)
)

combo_long <- combo_long[is.finite(combo_long$sec) & combo_long$sec > 0 &
                           is.finite(combo_long$n) & combo_long$n > 0, , drop = FALSE]

if (nrow(combo_long) == 0) {
  warning("No combined runtime data available. Skipping runtime_scaling_with_popkin plot.")
} else {
  xlim_combo <- range(combo_long$n, finite = TRUE)
  ylim_combo <- range(combo_long$sec, finite = TRUE)
  
  p_combo <- ggplot(combo_long, aes(x = n, y = sec, color = series)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_log10(limits = xlim_combo) +
    scale_y_log10(limits = ylim_combo) +
    scale_color_manual(values = c(
      "Popkin total (build Phi / A_min)" = COL_POPKIN["popkin"],
      "RSpectra (explicit Phi) total" = COL_RSPECTRA,
      "rARPACK (operator) total" = COL_RARPACK,
      "Explicit pipeline total = Popkin + RSpectra" = "black"
    )) +
    labs(
      title = "Runtime comparison (log10 axes): explicit Phi pipeline vs operator",
      x = "Sample size (n)",
      y = "Runtime (seconds, log10 axes)",
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
  
  ggsave(
    filename = "output/runtime_scaling_with_popkin.pdf",
    plot = p_combo,
    width = 7.6,
    height = 5.3
  )
} 

cat("Wrote PDFs to output/:\n",
    " - runtime_scaling_full.pdf\n",
    " - popkin_timing_breakdown.pdf\n",
    " - load_time_X_vs_Phi.pdf\n",
    " - runtime_scaling_with_popkin.pdf\n", sep = "")


