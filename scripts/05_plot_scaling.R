# 05_plot_scaling.R
# ------------------------------------------------------------
# Description: 
#   Reconstructs the full scaling curve from individual runs.  
#   
# ============================================================

files <- list.files("output", pattern="runtime_scaling_n_.*\\.csv", full.names=TRUE)
res <- do.call(rbind, lapply(files, read.csv))

pdf("output/runtime_scaling_full.pdf", width=7, height=5)
plot(res$n, res$rARPACK, type="b", pch=19, log="xy", col="red",
     xlab="Sample size (n)", ylab="Runtime (seconds)",
     main="Runtime Scaling: RSpectra vs rARPACK")
lines(res$n, res$RSpectra, type="b", pch=19, col="blue")
legend("topleft", legend=c("rARPACK","RSpectra"),
       col=c("red","blue"), pch=19)
dev.off()