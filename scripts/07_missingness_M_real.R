# scripts/07_missingness_M_real.R
# ------------------------------------------------------------
# Explore missingness structure in real PLINK genotype data
#
# Outputs:
#   - Histograms of missingness per locus / per individual
#   - Popkin missingness matrix M and approximations
#   - RMSD / relative-error summaries
#   - Scatterplots: true vs estimated pairwise overlap / missingness
#     with diagonal vs off-diagonal points distinguished
#   - Heatmaps of true missingness, approximations, residuals, kinship
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(genio)
  library(popkin)
  library(ggplot2)
  library(scales)
  library(grid)
})

# 0. Paths --------------------------------------------------------------------
dir <- "/Users/izzi/Desktop/Duke/Master's Project/Code/admix-popkin-masters/hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1"
plink_prefix <- file.path(dir, "hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1")
out_dir <- "/Users/izzi/Desktop/Duke/Master's Project/Code/admix-popkin-masters/output/missingness_hgdp"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- function(fname) file.path(out_dir, fname)

# 1. Load PLINK data ----------------------------------------------------------
cat("Reading PLINK prefix:\n", plink_prefix, "\n\n")

dat <- read_plink(plink_prefix)
X <- dat$X  # loci x individuals

m <- nrow(X)
n <- ncol(X)

cat("Loaded X with dimensions:", m, "loci x", n, "individuals\n")

# 2. Missingness diagnostics --------------------------------------------------
miss_locus <- rowMeans(is.na(X))
miss_indiv <- colMeans(is.na(X))

binwidth <- 0.0009

df_locus <- data.frame(missingness = miss_locus)
df_indiv <- data.frame(missingness = miss_indiv)

p1 <- ggplot(df_locus, aes(x = missingness)) +
  geom_histogram(
    binwidth = binwidth,
    fill = "#D55E00",
    color = "white",
    alpha = 0.9
  ) +
  labs(
    title = "Missingness per Locus",
    x = "Proportion of Individuals With Missing Genotype per Locus",
    y = "Count of Loci"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  out_path("hist_missingness_per_locus.png"),
  p1, width = 8, height = 5, dpi = 300
)

p2 <- ggplot(df_indiv, aes(x = missingness)) +
  geom_histogram(
    binwidth = binwidth,
    fill = "#0072B2",
    color = "white",
    alpha = 0.9
  ) +
  labs(
    title = "Missingness per Individual",
    x = "Proportion of Loci Missing per Individual",
    y = "Count of Individuals"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  out_path("hist_missingness_per_individual.png"),
  p2, width = 8, height = 5, dpi = 300
)

locus_summary <- data.frame(
  stat = names(summary(miss_locus)),
  value = as.numeric(summary(miss_locus))
)
indiv_summary <- data.frame(
  stat = names(summary(miss_indiv)),
  value = as.numeric(summary(miss_indiv))
)

write.csv(locus_summary, out_path("summary_missingness_locus.csv"), row.names = FALSE)
write.csv(indiv_summary, out_path("summary_missingness_individual.csv"), row.names = FALSE)

# 3. Popkin outputs including M -----------------------------------------------
cat("Running popkin(X, want_M = TRUE)...\n")
out <- popkin(X, want_M = TRUE)

K <- out$kinship
M <- out$M

# 4. Approximations to M ------------------------------------------------------
# Factorized approximation based on per-individual observed proportion
p <- colMeans(!is.na(X))
M_approx <- m * tcrossprod(p)

# Rank-1 eigen approximation
eig <- eigen(M, symmetric = TRUE)
lambda1 <- eig$values[1]
v1 <- eig$vectors[, 1]
M_eigen <- lambda1 * tcrossprod(v1)

# 5. RMSD and relative-error summaries ----------------------------------------
rmsd <- function(A, B) sqrt(mean((A - B)^2))

pct_rmsd <- function(A, B) {
  denom <- mean(A)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  100 * sqrt(mean((A - B)^2)) / denom
}

# Missingness scale
N_true   <- 1 - (M / m)
N_approx <- 1 - (M_approx / m)
N_eigen  <- 1 - (M_eigen / m)

metrics <- data.frame(
  comparison = c(
    "RMSD(M, M_approx)",
    "RMSD(M, M_eigen)",
    "RMSD(M_approx, M_eigen)",
    "RMSD(N_true, N_approx)",
    "RMSD(N_true, N_eigen)"
  ),
  value = c(
    rmsd(M, M_approx),
    rmsd(M, M_eigen),
    rmsd(M_approx, M_eigen),
    rmsd(N_true, N_approx),
    rmsd(N_true, N_eigen)
  )
)

metrics_pct <- data.frame(
  comparison = c(
    "Percent RMSD(M, M_approx)",
    "Percent RMSD(M, M_eigen)",
    "Percent RMSD(N_true, N_approx)",
    "Percent RMSD(N_true, N_eigen)"
  ),
  value = c(
    pct_rmsd(M, M_approx),
    pct_rmsd(M, M_eigen),
    pct_rmsd(N_true, N_approx),
    pct_rmsd(N_true, N_eigen)
  )
)

print(metrics)
print(metrics_pct)

write.csv(metrics, out_path("M_RMSD_metrics.csv"), row.names = FALSE)
write.csv(metrics_pct, out_path("M_relative_error_metrics.csv"), row.names = FALSE)

# 6. Pairwise scatterplots: true vs estimated ---------------------------------
idx <- upper.tri(M, diag = TRUE)

idx_rc <- which(idx, arr.ind = TRUE)
pair_type <- ifelse(idx_rc[,1] == idx_rc[,2], "Diagonal", "Off-diagonal")

df_pairs <- data.frame(
  true_M    = M[idx],
  approx_M  = M_approx[idx],
  eigen_M   = M_eigen[idx],
  pair_type = factor(pair_type, levels = c("Off-diagonal","Diagonal"))
)

df_pairs$true_N   <- 1 - (df_pairs$true_M / m)
df_pairs$approx_N <- 1 - (df_pairs$approx_M / m)
df_pairs$eigen_N  <- 1 - (df_pairs$eigen_M / m)

# Long format
df_M_long <- rbind(
  data.frame(true=df_pairs$true_M, estimate=df_pairs$approx_M,
             method="Factorized", pair_type=df_pairs$pair_type),
  data.frame(true=df_pairs$true_M, estimate=df_pairs$eigen_M,
             method="Eigen rank-1", pair_type=df_pairs$pair_type)
)

df_N_long <- rbind(
  data.frame(true=df_pairs$true_N, estimate=df_pairs$approx_N,
             method="Factorized", pair_type=df_pairs$pair_type),
  data.frame(true=df_pairs$true_N, estimate=df_pairs$eigen_N,
             method="Eigen rank-1", pair_type=df_pairs$pair_type)
)

# Randomize rows
set.seed(1)
df_M_long <- df_M_long[sample(nrow(df_M_long)), ]
df_N_long <- df_N_long[sample(nrow(df_N_long)), ]

# High contrast colors
approx_cols <- c(
  "Eigen rank-1" = "#B2182B",
  "Factorized"   = "#2166AC"
)

# ----- Scatterplot: M scale -----
p_M <- ggplot(df_M_long, aes(x=true, y=estimate)) +
  
  # background cloud
  geom_point(
    data=subset(df_M_long, pair_type=="Off-diagonal"),
    aes(color=method),
    shape=16,
    alpha=0.15,
    size=0.7
  ) +
  
  # diagonal triangles
  geom_point(
    data=subset(df_M_long, pair_type=="Diagonal"),
    aes(color=method),
    shape=17,
    size=1.5
  ) +
  
  geom_abline(slope=1, intercept=0, linetype=2) +
  
  scale_color_manual(
    values=approx_cols,
    name="Approximation"
  ) +
  
  labs(
    title="Pairwise Non-Missing Counts: True vs Estimated",
    subtitle="Diagonal entries shown as triangles",
    x="True M value",
    y="Estimated M value"
  ) +
  
  theme_bw(base_size=14) +
  theme(
    plot.title = element_text(face="bold"),
    plot.subtitle = element_text(size=12),
    legend.title = element_text(face="bold")
  )

ggsave(
  out_path("scatter_trueM_vs_estimates.png"),
  p_M, width=8.5, height=6, dpi=300
)

# ----- Scatterplot: missingness scale -----
p_N <- ggplot(df_N_long, aes(x=true, y=estimate)) +
  
  geom_point(
    data=subset(df_N_long, pair_type=="Off-diagonal"),
    aes(color=method),
    shape=16,
    alpha=0.15,
    size=0.7
  ) +
  
  geom_point(
    data=subset(df_N_long, pair_type=="Diagonal"),
    aes(color=method),
    shape=17,
    size=1.5
  ) +
  
  geom_abline(slope=1, intercept=0, linetype=2) +
  
  scale_color_manual(
    values=approx_cols,
    name="Approximation"
  ) +
  
  scale_x_continuous(labels = percent_format(accuracy=0.1)) +
  scale_y_continuous(labels = percent_format(accuracy=0.1)) +
  
  labs(
    title="Pairwise Missingness Proportion: True vs Estimated",
    subtitle="Diagonal entries shown as triangles",
    x="True missingness proportion",
    y="Estimated missingness proportion"
  ) +
  
  theme_bw(base_size=14) +
  theme(
    plot.title = element_text(face="bold"),
    plot.subtitle = element_text(size=12),
    legend.title = element_text(face="bold")
  )

ggsave(
  out_path("scatter_trueMissingness_vs_estimates.png"),
  p_N, width=8.5, height=6, dpi=300
)


# 7. Heatmaps -----------------------------------------------------------------
# Residuals in missingness scale
R_approx <- N_true - N_approx
R_eigen  <- N_true - N_eigen

# Sort using true M (mean pair overlap -> missingness score)
mean_pair_overlap <- rowMeans(M)
miss_score <- 1 - (mean_pair_overlap / m)
ord <- order(miss_score, decreasing = TRUE)

N_true_s   <- N_true[ord, ord]
N_approx_s <- N_approx[ord, ord]
N_eigen_s  <- N_eigen[ord, ord]
R_approx_s <- R_approx[ord, ord]
R_eigen_s  <- R_eigen[ord, ord]
K_s        <- K[ord, ord]

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
  
  if (!is.null(main) && !("main" %in% fmls)) {
    title(main, cex.main = 1.1)
  }
  
  box()
}

png(out_path("heatmap_N_true_sorted.png"), width = 1600, height = 1400, res = 200)
call_plot_popkin(
  N_true_s,
  main = "Missingness matrix: N_true = 1 - (M / m)\n(sorted by missingness from true M; higher = more missingness)",
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
  main = "Popkin kinship / coancestry (sorted by missingness order)",
  zlim = zlim_K
)
dev.off()

write.csv(
  data.frame(
    order_index = seq_along(ord),
    individual_index = ord,
    miss_score = miss_score[ord]
  ),
  out_path("ordering_by_trueM_missingness.csv"),
  row.names = FALSE
)

cat("Done. Outputs saved to:\n", out_dir, "\n")





