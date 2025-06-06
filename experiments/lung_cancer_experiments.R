library(diceRplus)
data(lung_cancer)


# Calcute internal_metrics based on consensus_evaluate of diceR
calculate_internal_metrics <- function(data,cluster_labels) {
  num_labels <- length(cluster_labels)

  # Create structure for consensus_evaluate
  cc_data <- array(cluster_labels, dim = c(num_labels, 1, 1, 1))
  # row_names <- rownames(data)
  row_names <- if (!is.null(rownames(data))) rownames(data) else seq_len(nrow(data))
    dimnames(cc_data) <- list(
    row_names,  # Primer nivel de nombres: nombres de las filas de Meat$x
    "R1",       # Repetition
    "RPGMMClu", # Clustering algorithm
    "5"         # Number of clusters
  )
  # Evaluation
  result_evaluation <- consensus_evaluate(data, cc_data)
  metrics_df <- result_evaluation$ii[[1]]
  metrics_list <- lapply(metrics_df[1, -1], function(x) unname(x))
  return(metrics_list)
}


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

scores <- list()
metrics_list <- list()

for (method in names(UFS_results$Results)) {
  print(method)
  top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]
  execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(lung_cancer$x[, top_features],
                                                                      lung_cancer$y,
                                                                      g=5,
                                                                      B=B,
                                                                      B.star=B.star,
                                                                      verb=TRUE))["elapsed"]
  # Calculate internal metrics
  internal_metrics <- calculate_internal_metrics(lung_cancer$x[, top_features],
                                                 out.clu_baseline_UFS$ensemble$label.vec)

  cat("ARI:", out.clu_baseline$ensemble$ari, "\n")
  cat("BIC:", out.clu_baseline$ensemble$bic.final, "\n")
  print(internal_metrics)
  metrics_list[[method]] <- internal_metrics
}

# Convertir cada columna a numérica
metrics_df <- as.data.frame(lapply(metrics_df, function(x) as.numeric(unlist(x))))

# Restaurar nombres de los métodos
rownames(metrics_df) <- names(metrics_list)

# Ahora sí, aplicar la normalización
normalize <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
metrics_norm <- as.data.frame(lapply(metrics_df, normalize))

# Restaurar nombres de los métodos en la tabla normalizada
rownames(metrics_norm) <- rownames(metrics_df)

metrics_norm[is.na(metrics_norm)] <- 0

# Invertir métricas donde menor es mejor
metrics_norm$davies_bouldin <- 1 - metrics_norm$davies_bouldin
metrics_norm$Compactness <- 1 - metrics_norm$Compactness
metrics_norm$Connectivity <- 1 - metrics_norm$Connectivity
metrics_norm$c_index <- 1 - metrics_norm$c_index
metrics_norm$mcclain_rao <- 1 - metrics_norm$mcclain_rao
metrics_norm$sd_dis <- 1 - metrics_norm$sd_dis
metrics_norm$ray_turi <- 1 - metrics_norm$ray_turi



# Calcular score después de la normalización
weights <- c(
  calinski_harabasz = 0.2,
  dunn = 0.15,
  pbm = 0.2,
  tau = 0,
  gamma = 0,
  c_index = 0,
  mcclain_rao = 0,
  sd_dis = 0,
  ray_turi = 0,
  g_plus = 0,
  silhouette = 0.15,
  s_dbw = 0,
  davies_bouldin = 0.1,
  Compactness = 0.1,
  Connectivity = 0.1
)

scores <- rowSums(sweep(metrics_norm, 2, weights, "*"))  # Multiplica por pesos

# Obtener el mejor método según el score
best_method <- names(which.max(scores))

cat("El mejor método basado en métricas internas es:", best_method, "\n")


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
  internal_metrics <- calculate_internal_metrics(lung_cancer$x_filtered,
                                                     out.clu_baseline_UFS$ensemble$label.vec)

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
  print(internal_metrics)

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
