library(diceRplus)
data(lung_cancer)


#lung_cancer$x <- scale(lung_cancer$x)

# Verify normalizacion
mean(lung_cancer$x)
sd(lung_cancer$x)


cor_mat <- cor(lung_cancer$x, use = "pairwise.complete.obs")
mean_cor <- mean(abs(cor_mat[upper.tri(cor_mat)]))

var_vec <- apply(lung_cancer$x, 2, var, na.rm = TRUE)
low_var_prop <- mean(var_vec < quantile(var_vec, 0.1))

pca <- prcomp(lung_cancer$x, scale. = TRUE)
var_exp <- summary(pca)$importance[2, ]
cum_var_90 <- which(cumsum(var_exp) > 0.9)[1]

ufs_candidates <- c()

if(mean_cor > 0.5) ufs_candidates <- c(ufs_candidates, "MCFS", "InfFS", "Laplacian")
if(low_var_prop > 0.3) ufs_candidates <- c(ufs_candidates, "UDFS", "NDFS", "InfFS")
if(cum_var_90 < ncol(lymphoma$x) * 0.2) ufs_candidates <- c(ufs_candidates, "UDFS", "NDFS", "InfFS")
if(ncol(lung_cancer$x) > 3000) ufs_candidates <- c(ufs_candidates, "Laplacian", "SPEC", "InfFS")


# Cover at least 3 types of UFS
ufs_candidates <- unique(c(ufs_candidates, "UDFS", "InfFS", "Laplacian", "NDFS"))
ufs_candidates <- unique(ufs_candidates)
print(ufs_candidates)
# Parameter Num clusters
# UDFS: Yes
# NDFS: No
# InfFS: No
# Laplacian: No
# SPEC: No


UFS_Methods <- list(
  "InfFS" = FALSE,
  "Laplacian" = TRUE,
  # "MCFS" = TRUE,
  "LLCFS" = TRUE,
  "CFS" = FALSE,
  "FSASL" = TRUE,
  "DGUFS" = TRUE,
  "UFSOL" = TRUE,
  "SPEC" = TRUE,
  # "SOCFS" = TRUE, # Won't work
  # "SOGFS" = TRUE,
  "UDFS" = TRUE,
  "SRCFS" = TRUE,
  # "FMIUFS" = TRUE,
  # "UAR_HKCMI" = TRUE,
  "RNE" = TRUE,
  # "FRUAFR" = TRUE, # Won't work
  # "U2FS" = TRUE,
  "RUFS" = TRUE,
  "NDFS" = TRUE,
  "EGCFS" = TRUE,
  "CNAFS" = TRUE,
  "Inf-FS2020" = TRUE
)

# Run just once
#UFS_results <- runUFS_manual_params(lung_cancer$x, ufs_candidates)
#save(UFS_results, file = "experiments/UFS_lung_cancer.RData")


load("experiments/UFS_lung_cancer.RData")

N <- cum_var_90
B <- N
B.star <- round(B/10)

cat("N:", N, "\nB:", B, "\nB.star:", B.star, "\n")


for (method in names(UFS_results$Results)) {
  print(method)
  top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]
  execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(lung_cancer$x[, top_features],
                                                                      lung_cancer$y,
                                                                      g=5,
                                                                      B=B,
                                                                      B.star=B.star,
                                                                      verb=TRUE))["elapsed"]
  print(out.clu_baseline$ensemble)
}

# Fixed UFS: Lapacian
method <- names(UFS_results$Results)[4]
print(method)

# B <- 500
# B.star <- 50


# Iterate over different values of N (from 20 to 100 in steps of 10)
for (N in seq(20, cum_var_90*4, by = 10)) {
  print(N)
  top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]

  # Run RPGMMClu_parallel
  execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(lung_cancer$x[, top_features],
                                                                          lung_cancer$y,
                                                                          g = 5,
                                                                          B = B,
                                                                          B.star = B.star))["elapsed"]
  # Calculate internal metrics
  #internal_metrics_UFS <- calculate_internal_metrics(lung_cancer$x_filtered,
  #                                                   out.clu_baseline_UFS$ensemble$label.vec)

  # Log the experiment
  exp_data <- experiment_logger(
    description = paste("Clustering with RPGMMClu - UFS:", method, "- Features:", N),
    dataset = "Meat",
    ensemble_method = "RPGMMClu_parallel",
    ensemble_method_params = list(g = 5, B = B, B.star = B.star),
    UFS_method = method,
    UFS_method_params = list(default = 'YES'),
    num_features = N,
    features = top_features,
    execution_time = as.numeric(execution_time),
    labels_clustering = out.clu_baseline_UFS$ensemble$label.vec,
    #internal_metrics = internal_metrics_UFS,
    external_metrics = list(ensemble_ari = out.clu_baseline_UFS$ensemble$ari[[1]])
  )
  cat("N:", N, "\nARI:", out.clu_baseline$ensemble$ari, "\n")
  # save_experiment(exp_data)
}



# Fixed UFS: InfFS
method <- names(UFS_results$Results)[1]
print(method)
N=1000
top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]

B <- 500
B.star=50

# Run RPGMMClu
execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(lung_cancer$x[, top_features],
                                                                    lung_cancer$y,
                                                                    g=5,
                                                                    B=B,
                                                                    B.star=B.star,
                                                                    verb=TRUE))["elapsed"]
