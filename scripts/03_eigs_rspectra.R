# 03_eigs_rspectra.R

suppressPackageStartupMessages({
  library(argparser)
  library(RSpectra)
})

p <- arg_parser("Top-k eigendecomposition via RSpectra (explicit matrix)")
p <- add_argument(p, "--phi", help="Path to Phi/Theta matrix tsv", default="output/Phi.tsv")
p <- add_argument(p, "--k", help="Number of eigenvectors", type="integer", default=10)
p <- add_argument(p, "--out", help="Output prefix", default="output/eigs_rspectra")
argv <- parse_args(p)

Phi <- as.matrix(read.table(argv$phi, header=TRUE, row.names=1))

t0 <- Sys.time()
fit <- eigs_sym(Phi, k=argv$k, which="LA")
runtime <- as.numeric(difftime(Sys.time(), t0, units="secs"))

dir.create(dirname(argv$out), recursive=TRUE, showWarnings=FALSE)

write.table(fit$values, paste0(argv$out, "_values.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(fit$vectors, paste0(argv$out, "_vectors.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# runtime 
write.csv(data.frame(runtime_sec = runtime),
          paste0(argv$out, "_runtime.csv"),
          row.names = FALSE)

cat("RSpectra runtime (sec):", runtime, "\n")

