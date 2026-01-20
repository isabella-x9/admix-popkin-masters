# scripts/04c_run_rspectra.R
suppressPackageStartupMessages({
  library(RSpectra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/04c_run_rspectra.R <n>")

n <- as.integer(args[1])
k <- min(10, n - 1)

phi_path   <- sprintf("output/tmp/Phi_n_%d.tsv", n)
out_prefix <- sprintf("output/tmp/rspectra_n_%d", n)

rt_path    <- paste0(out_prefix, "_runtime.csv")
time_path  <- paste0(out_prefix, "_timing.csv")

stopifnot(file.exists(phi_path))
dir.create("output/tmp", showWarnings = FALSE, recursive = TRUE)

# Skip if runtime exists 
if (file.exists(rt_path) && file.exists(time_path)) {
  cat("RSpectra outputs exist, skipping:", rt_path, "\n")
  quit(save = "no", status = 0)
}

t_total0 <- Sys.time()

# Load Phi
t0 <- Sys.time()
Phi <- as.matrix(read.table(phi_path, header = TRUE, row.names = 1))
t_loadPhi <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Eigs
t0 <- Sys.time()
fit <- eigs_sym(Phi, k = k, which = "LA")
t_eigs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Write outputs
t0 <- Sys.time()
write.table(fit$values, paste0(out_prefix, "_values.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(fit$vectors, paste0(out_prefix, "_vectors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
t_write <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t_total <- as.numeric(difftime(Sys.time(), t_total0, units = "secs"))

# Save runtime (keep existing file)
write.csv(data.frame(runtime_sec = t_eigs), rt_path, row.names = FALSE)

# Save detailed timing
write.csv(
  data.frame(
    n = n,
    k = k,
    load_Phi_sec = t_loadPhi,
    eigs_sec = t_eigs,
    write_outputs_sec = t_write,
    total_sec = t_total
  ),
  time_path,
  row.names = FALSE
)

cat("RSpectra runtime (sec):", t_eigs, "\n")
cat("RSpectra done.\n")
cat("Wrote:", time_path, "\n")

