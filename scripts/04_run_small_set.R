ns <- c(10, 100, 1000)

for (n in ns) {
  cat("\n=== n =", n, "===\n")
  system(sprintf("Rscript scripts/04a_simulate_geno.R %d", n))
  system(sprintf("Rscript scripts/04b_compute_phi_amin.R %d", n))
  system(sprintf("Rscript scripts/04c_run_rspectra.R %d", n))
  system(sprintf("Rscript scripts/04d_run_rarpack.R %d", n))
  system(sprintf("Rscript scripts/04e_combine_csv.R %d", n))
}
