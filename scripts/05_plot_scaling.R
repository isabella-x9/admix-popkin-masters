# 05_plot_scaling.R
# ------------------------------------------------------------
# Description:
#   Plot scaling for small demo runs only: n = 10, 100, 1000
#   Reads Stage 04 outputs: output/runtime_scaling_n_<n>.csv
# ============================================================

target_n <- c(10, 100, 1000)

# load exactly the three per-n files
files <- file.path("output", paste0("runtime_scaling_n_", target_n, ".csv"))
missing <- files[!file.exists(files)]
if (length(missing)) stop("Missing files:\n", paste(missing, collapse = "\n"))

res <- do.call(rbind, lapply(files, read.csv))
res <- res[match(target_n, res$n), ]

# save combined
write.csv(res, "output/runtime_scaling_small_combined.csv", row.names = FALSE)

# plot (log–log), allow RSpectra NA
ok_ra <- is.finite(res$rARPACK)  & res$rARPACK  > 0
ok_rs <- is.finite(res$RSpectra) & res$RSpectra > 0
yr <- range(c(res$rARPACK[ok_ra], res$RSpectra[ok_rs]), finite = TRUE)

pdf("output/runtime_scaling_small.pdf", width = 6, height = 5)

plot(res$n[ok_ra], res$rARPACK[ok_ra],
     type = "b", pch = 19, log = "xy",
     col = "firebrick", lwd = 2,
     xlab = "Sample size (n)",
     ylab = "Runtime (seconds)",
     main = "Runtime Scaling (n = 10, 100, 1000)",
     ylim = yr)

if (any(ok_rs)) {
  lines(res$n[ok_rs], res$RSpectra[ok_rs],
        type = "b", pch = 19,
        col = "steelblue", lwd = 2)
}

legend("topleft",
       legend = c("rARPACK", "RSpectra"),
       col = c("firebrick", "steelblue"),
       pch = 19,
       lwd = 2)

dev.off()

cat("Saved:\n",
    " - output/runtime_scaling_small.pdf\n",
    " - output/runtime_scaling_small_combined.csv\n", sep = "")
