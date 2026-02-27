# scripts/07_missingness_M_real.R 
# ------------------------------------------------------------ 
# Explore missingness structure in real PLINK genotype data

library(genio)
library(popkin)
library(ggplot2)

# 0. Load paths ----------------------------------------------------------------
dir <- "/Users/izzi/Desktop/Duke/Master's Project/Code/admix-popkin-masters/hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1"
plink_prefix <- file.path(dir, "hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1")
out_dir <- "/Users/izzi/Desktop/Duke/Master's Project/Code/admix-popkin-masters/output/missingness_hgdp"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- function(fname) file.path(out_dir, fname)


# 1. Load PLINK data -----------------------------------------------------------
cat("Reading PLINK prefix:\n", plink_prefix, "\n\n")

dat <- read_plink(plink_prefix)
X <- dat$X  # loci x individuals

m <- nrow(X)
n <- ncol(X)

cat("Loaded X with dimensions:", m, "loci x", n, "individuals\n")


# 2. Missingness diagnostics ---------------------------------------------------
miss_locus <- rowMeans(is.na(X))
miss_indiv <- colMeans(is.na(X))

## Histograms with fixed binwidth
binwidth <- 0.0009

df_locus <- data.frame(missingness = miss_locus)
df_indiv <- data.frame(missingness = miss_indiv)

p1 <- ggplot(df_locus, aes(x = missingness)) +
  geom_histogram(binwidth = binwidth,
                 fill = "#D55E00",
                 color = "white",
                 alpha = 0.9) +
  labs(
    title = "Missingness per Locus",
    subtitle = "Proportion of individuals with missing genotype per locus",
    x = "Proportion Missing",
    y = "Count of Loci"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(out_path("hist_missingness_per_locus.png"),
       p1, width = 8, height = 5, dpi = 300)

p2 <- ggplot(df_indiv, aes(x = missingness)) +
  geom_histogram(binwidth = binwidth,
                 fill = "#0072B2",
                 color = "white",
                 alpha = 0.9) +
  labs(
    title = "Missingness per Individual",
    subtitle = "Proportion of loci missing per individual",
    x = "Proportion Missing",
    y = "Count of Individuals"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(out_path("hist_missingness_per_individual.png"),
       p2, width = 8, height = 5, dpi = 300)

## Save missingness summaries (CSV)
locus_summary <- data.frame(
  stat = names(summary(miss_locus)),
  value = as.numeric(summary(miss_locus))
)
indiv_summary <- data.frame(
  stat = names(summary(miss_indiv)),
  value = as.numeric(summary(miss_indiv))
)

write.csv(locus_summary, 
          out_path("summary_missingness_locus.csv"), 
          row.names = FALSE)
write.csv(indiv_summary, 
          out_path("summary_missingness_individual.csv"), 
          row.names = FALSE)


# 3. Popkin outputs including M ------------------------------------------------
cat("Running popkin(X, want_M = TRUE)...\n")
out <- popkin(X, want_M = TRUE) 

## kinship estimate (for reference)
K <- out$kinship 

## matrix of non-missing pair counts
M <- out$M 


# 4. Factorized approximation M_approx = m * tcrossprod(p) ---------------------
p <- colMeans(!is.na(X))  # proportion observed per person
M_approx <- m * tcrossprod(p)


# 5. Theoretically optimal rank-1 approximation via eigen(M) -------------------
eig <- eigen(M, symmetric = TRUE)
lambda1 <- eig$values[1]
v1 <- eig$vectors[, 1]
M_eigen <- lambda1 * tcrossprod(v1)


# 6. Compare using RMSD
rmsd <- function(A, B) sqrt(mean((A - B)^2))

metrics <- data.frame(
  comparison = c("RMSD(M, M_approx)", "RMSD(M, M_eigen)", "RMSD(M_approx, M_eigen)"),
  value = c(
    rmsd(M, M_approx),
    rmsd(M, M_eigen),
    rmsd(M_approx, M_eigen)
  )
)

print(metrics)
write.csv(metrics, out_path("M_RMSD_metrics.csv"), row.names = FALSE)


# 7. Visualize as heatmaps (multi-panel)

# Normalize to missingness fraction: higher = more missingness
N_true   <- 1 - (M / m)
N_approx <- 1 - (M_approx / m)
N_eigen  <- 1 - (M_eigen / m)

# Residuals in missingness scale
R_approx <- N_true - N_approx
R_eigen <- N_true - N_eigen

# Sort using TRUE M (mean pair overlap -> missingness score)
mean_pair_overlap <- rowMeans(M)
miss_score <- 1 - (mean_pair_overlap / m)
ord <- order(miss_score, decreasing = TRUE)

# Apply ordering consistently
N_true_s <- N_true[ord, ord]
N_approx_s <- N_approx[ord, ord]
N_eigen_s <- N_eigen[ord, ord]
R_approx_s <- R_approx[ord, ord]
R_eigen_s <- R_eigen[ord, ord]
K_s <- K[ord, ord]

# Consistent scales
zlim_N <- range(c(N_true_s, N_approx_s, N_eigen_s), finite = TRUE)

lim_R  <- max(abs(c(R_approx_s, R_eigen_s)), finite = TRUE)
zlim_R <- c(-lim_R, lim_R)

zlim_K <- range(K_s, finite = TRUE)


call_plot_popkin <- function(A, main = NULL, zlim = NULL) {
  fmls <- names(formals(plot_popkin))
  args <- list(A)
  
  if (!is.null(main) && "main" %in% fmls) args$main <- main
  if (!is.null(zlim) && "zlim" %in% fmls) args$zlim <- zlim
  
  do.call(plot_popkin, args)
  
  # If plot_popkin doesn't accept main, add it after
  if (!is.null(main) && !("main" %in% fmls)) {
    title(main, cex.main = 1.1)
  }
  
  # Make axes/titles a bit cleaner
  box()
}

# Save plots
png(out_path("heatmap_N_true_sorted.png"), width = 1600, height = 1400, res = 200)
call_plot_popkin(
  N_true_s,
  main = "Missingness matrix: N_true = 1 - (M/m)\n(sorted by missingness from true M; higher = more missingness)",
  zlim = zlim_N
)
dev.off()

png(out_path("heatmap_N_approx_sorted.png"), width = 1600, height = 1400, res = 200)
call_plot_popkin(
  N_approx_s,
  main = "Missingness approximation: N_approx from M_approx = m * p p^T\n(same ordering; higher = more missingness)",
  zlim = zlim_N
)
dev.off()

png(out_path("heatmap_N_eigen_sorted.png"), width = 1600, height = 1400, res = 200)
call_plot_popkin(
  N_eigen_s,
  main = "Missingness approximation: N_eigen from rank-1 eigen(M)\n(same ordering; higher = more missingness)",
  zlim = zlim_N
)
dev.off()

png(out_path("heatmap_residual_N_minus_Napprox_sorted.png"), width = 1600, height = 1400, res = 200)
call_plot_popkin(
  R_approx_s,
  main = "Residual (missingness scale): N_true - N_approx\n(sorted; symmetric scale)",
  zlim = zlim_R
)
dev.off()

png(out_path("heatmap_residual_N_minus_Neigen_sorted.png"), width = 1600, height = 1400, res = 200)
call_plot_popkin(
  R_eigen_s,
  main = "Residual (missingness scale): N_true - N_eigen\n(sorted; symmetric scale)",
  zlim = zlim_R
)
dev.off()

png(out_path("heatmap_kinship_sorted.png"), width = 1600, height = 1400, res = 200)
call_plot_popkin(
  K_s,
  main = "Popkin kinship/coancestry (sorted by missingness order)",
  zlim = zlim_K
)
dev.off()

# Save ordering used for all plots
write.csv(
  data.frame(order_index = seq_along(ord),
             individual_index = ord,
             miss_score = miss_score[ord]),
  out_path("ordering_by_trueM_missingness.csv"),
  row.names = FALSE
)

cat("Done. Outputs saved to:\n", out_dir, "\n")






