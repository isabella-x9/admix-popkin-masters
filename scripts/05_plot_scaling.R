# 05_plot_scaling.R
# ------------------------------------------------------------
# Plots:
#   1) Runtime scaling: RSpectra vs rARPACK 
#   2) Popkin step breakdown (load X / popkin / popkin_A / write)
#   3) Load-time comparison: load X vs load Phi
#
# Reads:
#   output/runtime_scaling_n_<n>.csv  (created by 04e)
#   output/tmp/popkin_n_<n>_timing.csv
#   output/tmp/rspectra_n_<n>_timing.csv
#   output/tmp/rarpack_n_<n>_timing.csv
# ------------------------------------------------------------

N_VALUES <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)

# ---------- helper ----------
safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path)
}

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ---------- combined totals table (forces all N_VALUES on x-axis) ----------
rows <- lapply(N_VALUES, function(n) {
  f <- sprintf("output/runtime_scaling_n_%d.csv", n)
  if (!file.exists(f)) {
    return(data.frame(n = n, m = NA_real_, RSpectraTotal = NA_real_, rARPACKTotal = NA_real_))
  }
  d <- read.csv(f)
  data.frame(
    n = d$n[1],
    m = d$m[1],
    RSpectraTotal = d$RSpectraTotal[1],
    rARPACKTotal  = d$rARPACKTotal[1]
  )
})
res <- do.call(rbind, rows)

write.csv(res, "output/runtime_scaling_combined.csv", row.names = FALSE)

ok_ra <- is.finite(res$rARPACKTotal)  & res$rARPACKTotal  > 0
ok_rs <- is.finite(res$RSpectraTotal) & res$RSpectraTotal > 0
yr <- range(c(res$rARPACKTotal[ok_ra], res$RSpectraTotal[ok_rs]), finite = TRUE)

pdf("output/runtime_scaling_full.pdf", width = 7, height = 5)

plot(res$n[ok_ra], res$rARPACKTotal[ok_ra],
     type = "b", pch = 19, log = "xy",
     xlab = "Sample size (n)",
     ylab = "Total runtime (seconds)",
     main = "Runtime Scaling: explicit (Phi) vs operator (no Phi)",
     ylim = yr)

if (any(ok_rs)) {
  lines(res$n[ok_rs], res$RSpectraTotal[ok_rs],
        type = "b", pch = 19)
}

legend("topleft",
       legend = c("rARPACK total", "RSpectra total"),
       pch = 19, lwd = 2)

dev.off()

# ---------- Popkin step timing breakdown (Phi-version bottleneck) ----------
pop_rows <- lapply(N_VALUES, function(n) {
  f <- sprintf("output/tmp/popkin_n_%d_timing.csv", n)
  d <- safe_read(f)
  if (is.null(d)) {
    return(data.frame(
      n = n,
      load_X_sec = NA_real_,
      popkin_sec = NA_real_,
      popkin_A_sec = NA_real_,
      write_outputs_sec = NA_real_,
      total_sec = NA_real_
    ))
  }
  d[, c("n","load_X_sec","popkin_sec","popkin_A_sec","write_outputs_sec","total_sec")]
})
pop <- do.call(rbind, pop_rows)

pdf("output/popkin_timing_breakdown.pdf", width = 7.5, height = 5)

ok <- is.finite(pop$total_sec) & pop$total_sec > 0
plot(pop$n[ok], pop$load_X_sec[ok],
     type = "b", pch = 19, log = "xy",
     xlab = "Sample size (n)",
     ylab = "Seconds (log-log)",
     main = "Popkin step breakdown (explicit-Phi pipeline)")

lines(pop$n[ok], pop$popkin_sec[ok], type = "b", pch = 19)
lines(pop$n[ok], pop$popkin_A_sec[ok], type = "b", pch = 19)
lines(pop$n[ok], pop$write_outputs_sec[ok], type = "b", pch = 19)

legend("topleft",
       legend = c("load X", "popkin()", "popkin_A() (A_min)", "write Phi/A_min"),
       pch = 19, lwd = 2)

dev.off()

# ---------- Load-time comparison: load X vs load Phi ----------
rs_rows <- lapply(N_VALUES, function(n) {
  f <- sprintf("output/tmp/rspectra_n_%d_timing.csv", n)
  d <- safe_read(f)
  if (is.null(d)) return(data.frame(n = n, load_Phi_sec = NA_real_))
  data.frame(n = n, load_Phi_sec = d$load_Phi_sec[1])
})
rs <- do.call(rbind, rs_rows)

ra_rows <- lapply(N_VALUES, function(n) {
  f <- sprintf("output/tmp/rarpack_n_%d_timing.csv", n)
  d <- safe_read(f)
  if (is.null(d)) return(data.frame(n = n, load_X_sec = NA_real_))
  data.frame(n = n, load_X_sec = d$load_X_sec[1])
})
ra <- do.call(rbind, ra_rows)

pdf("output/load_time_X_vs_Phi.pdf", width = 7.5, height = 5)

okx <- is.finite(ra$load_X_sec) & ra$load_X_sec > 0
okp <- is.finite(rs$load_Phi_sec) & rs$load_Phi_sec > 0
yr2 <- range(c(ra$load_X_sec[okx], rs$load_Phi_sec[okp]), finite = TRUE)

plot(ra$n[okx], ra$load_X_sec[okx],
     type = "b", pch = 19, log = "xy",
     xlab = "Sample size (n)",
     ylab = "Seconds",
     main = "I/O cost grows: loading X vs loading Phi",
     ylim = yr2)

if (any(okp)) {
  lines(rs$n[okp], rs$load_Phi_sec[okp], type = "b", pch = 19)
}

legend("topleft",
       legend = c("load X (rARPACK)", "load Phi (RSpectra)"),
       pch = 19, lwd = 2)

dev.off()

cat("Saved:\n",
    " - output/runtime_scaling_full.pdf\n",
    " - output/runtime_scaling_combined.csv\n",
    " - output/popkin_timing_breakdown.pdf\n",
    " - output/load_time_X_vs_Phi.pdf\n", sep = "")

