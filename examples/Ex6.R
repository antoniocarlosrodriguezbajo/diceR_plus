# 1. Función SNR multiclass (one-vs-all)
# 1. Carga y preparación de datos
data(leukemia)
dat <- leukemia[, !names(leukemia) %in% "class"]
ref.cl <- as.integer(leukemia$class)

# 2. Selección de genes por SNR y permutación
get_snr_and_pval <- function(mat, labels, nperm = 500) {
  classes <- unique(labels)
  snr_obs <- matrix(NA, nrow = length(classes), ncol = ncol(mat))
  pvals <- matrix(NA, nrow = length(classes), ncol = ncol(mat))
  colnames(snr_obs) <- colnames(mat)
  colnames(pvals) <- colnames(mat)
  rownames(snr_obs) <- rownames(pvals) <- classes
  for (k in seq_along(classes)) {
    cl <- classes[k]
    bin_labels <- ifelse(labels == cl, cl, paste0("not_", cl))
    idx1 <- which(bin_labels == cl)
    idx2 <- which(bin_labels != cl)
    mu1 <- colMeans(mat[idx1, , drop = FALSE])
    mu2 <- colMeans(mat[idx2, , drop = FALSE])
    sd1 <- apply(mat[idx1, , drop = FALSE], 2, sd)
    sd2 <- apply(mat[idx2, , drop = FALSE], 2, sd)
    snr0 <- (mu1 - mu2) / (sd1 + sd2)
    snr_obs[k, ] <- snr0

    snr_perm <- replicate(nperm, {
      perm <- sample(bin_labels)
      idx1p <- which(perm == cl)
      idx2p <- which(perm != cl)
      mu1p <- colMeans(mat[idx1p, , drop = FALSE])
      mu2p <- colMeans(mat[idx2p, , drop = FALSE])
      sd1p <- apply(mat[idx1p, , drop = FALSE], 2, sd)
      sd2p <- apply(mat[idx2p, , drop = FALSE], 2, sd)
      (mu1p - mu2p) / (sd1p + sd2p)
    })
    pvals[k, ] <- rowMeans(abs(snr_perm) >= abs(snr0))
  }
  list(snr = snr_obs, pval = pvals)
}
res <- get_snr_and_pval(as.matrix(dat), ref.cl, nperm = 500)
genes_sig <- colnames(dat)[colSums(res$pval < 0.05) > 0]
dat_selected <- dat[, genes_sig, drop = FALSE]

# 3. Normalización
normalize_rows_and_cols <- function(M) {
  M1 <- t(scale(t(M), center = TRUE, scale = TRUE)) # filas
  M2 <- scale(M1, center = TRUE, scale = TRUE)      # columnas
  return(M2)
}
dat_selected <- normalize_rows_and_cols(dat_selected)

# 4. Clustering consenso
library(diceR)
set.seed(123)
resultado <- dice(
  data      = dat_selected,
  nk        = 3:5,
  reps      = 500,
  algorithms= "hc",
  hc.method = "average",
  cons.funs = c("majority"),
  evaluate  = TRUE,
  ref.cl    = ref.cl,
  progress  = TRUE
)

# 5. Matriz de consenso y dendrograma (elige K según consenso)
cm <- consensus_matrix(resultado$Ecomp)
hc <- hclust(as.dist(1 - cm), method = "average")
cl_final <- cutree(hc, k = 3)

# 6. Evaluación
library(mclust)
adjustedRandIndex(cl_final, ref.cl)
library(aricode)
NMI(cl_final, ref.cl)

# 7. Visualización si quieres
graph_cdf(resultado)
graph_heatmap(resultado)
graph_delta_area(resultado)
