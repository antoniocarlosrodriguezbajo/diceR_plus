# Anderlucci's meat with UFS integrated into the package and internal metrics
# Anderlucci's: https://rdrr.io/cran/RPEClust/f/
# UFS: https://github.com/farhadabedinzadeh/AutoUFSTool

# Calcute internal_metrics based on consensus_evaluate of diceR
calculate_internal_metrics <- function(data,cluster_labels) {
  num_labels <- length(cluster_labels)

  # Create structure for consensus_evaluate
  cc_data <- array(cluster_labels, dim = c(num_labels, 1, 1, 1))
  row_names <- rownames(data)
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

# Meat dataset
data(Meat)

############################################################
# Experiment in Anderlucci's paper (original implementation)
############################################################
B <- 1000
B.star=100

# Run RPGMMClu
execution_time <- system.time(out.clu_baseline <- RPGMMClu(Meat$x,
                                                           Meat$y,
                                                           g=5,
                                                           B=B,
                                                           B.star=B.star,
                                                           verb=TRUE))["elapsed"]

# Calculate internal metrics
internal_metrics <- calculate_internal_metrics(Meat$x,
                                               out.clu_baseline$ensemble$label.vec)
# Log the experiment
exp_data <- experiment_logger(
  description = "Clustering with RPGMMClu - Experiment in Anderlucci's paper",
  dataset = "Meat",
  ensemble_method = "RPGMMClu",
  ensemble_method_params = list(g = 5, B = B, B.star = B.star),
  execution_time = as.numeric(execution_time),
  labels_clustering = out.clu_baseline$ensemble$label.vec,
  internal_metrics = internal_metrics,
  external_metrics = list(ensemble_ari = out.clu_baseline$ensemble$ari[[1]])
)
save_experiment(exp_data)

############################################################
# Experiment in Anderlucci's paper (parallel implementation)
############################################################
B <- 1000
B.star=100

# Run RPGMMClu_parallel
execution_time <- system.time(out.clu_baseline_p <- RPGMMClu_parallel(Meat$x,
                                                                      Meat$y,
                                                                      g=5,
                                                                      B=B,
                                                                      B.star=B.star))["elapsed"]

# Calculate internal metrics
internal_metrics_p <- calculate_internal_metrics(Meat$x,
                                               out.clu_baseline_p$ensemble$label.vec)
# Log the experiment
exp_data <- experiment_logger(
  description = "Clustering with RPGMMClu - Experiment in Anderlucci's paper (parallel implementation)",
  dataset = "Meat",
  ensemble_method = "RPGMMClu_parallel",
  ensemble_method_params = list(g = 5, B = B, B.star = B.star),
  execution_time = as.numeric(execution_time),
  labels_clustering = out.clu_baseline_p$ensemble$label.vec,
  internal_metrics = internal_metrics_p,
  external_metrics = list(ensemble_ari = out.clu_baseline_p$ensemble$ari[[1]])
)
save_experiment(exp_data)


############################################################
# Experiment with UFS
############################################################
UFS_Methods <- list(
  "InfFS" = FALSE,
  "Laplacian" = TRUE,
  "MCFS" = TRUE,
  "LLCFS" = TRUE,
  "CFS" = FALSE,
  "FSASL" = TRUE,
  "DGUFS" = TRUE,
  "UFSOL" = TRUE,
  "SPEC" = TRUE,
  # "SOCFS" = TRUE, # Won't work
  "SOGFS" = TRUE,
  "UDFS" = TRUE,
  "SRCFS" = TRUE,
  "FMIUFS" = TRUE,
  "UAR_HKCMI" = TRUE,
  "RNE" = TRUE,
  # "FRUAFR" = TRUE, # Won't work
  "U2FS" = TRUE,
  "RUFS" = TRUE,
  "NDFS" = TRUE,
  "EGCFS" = TRUE,
  "CNAFS" = TRUE,
  "Inf-FS2020" = TRUE
)
# Run just once
# UFS_Results <- runUFS(Meat$x, UFS_Methods)
# save(UFS_Results, file = "experiments/UFS_Meat.RData")

load("experiments/UFS_Meat.RData")

# Experiment with the 100 most relevant features

N <- 100
UFS_Results$Results$"Inf-FS2020"$Result[[1]][1:N]
UFS_Results$ExecutionTimes$"Inf-FS2020"

for (method in names(UFS_Results$Results)) {
  print(method)
  top_features <- UFS_Results$Results[[method]]$Result[[1]][1:N]
  Meat$x_filtered <- Meat$x[, top_features]
  # Run RPGMMClu_parallel
  execution_time <- system.time(out.clu_baseline_UFS <- RPGMMClu_parallel(Meat$x_filtered,
                                                                          Meat$y,
                                                                          g=5,
                                                                          B=B,
                                                                          B.star=B.star))["elapsed"]
  # Calculate internal metrics
  internal_metrics_UFS <- calculate_internal_metrics(Meat$x,
                                                     out.clu_baseline_UFS$ensemble$label.vec)
  # Log the experiment
  exp_data <- experiment_logger(
    description = "Clustering with RPGMMClu - UFS:Inf-FS2020",
    dataset = "Meat",
    ensemble_method = "RPGMMClu_parallel",
    ensemble_method_params = list(g = 5, B = B, B.star = B.star),
    execution_time = as.numeric(execution_time),
    labels_clustering = out.clu_baseline_UFS$ensemble$label.vec,
    internal_metrics = internal_metrics_UFS,
    external_metrics = list(ensemble_ari = out.clu_baseline_UFS$ensemble$ari[[1]])
  )
  save_experiment(exp_data)
}

############################

top_features <- UFS_Results$Results$"Inf-FS2020"$Result[[1]][1:N]
Meat$x_filtered <- Meat$x[, top_features]

# Run RPGMMClu_parallel
execution_time <- system.time(out.clu_baseline_UFS <- RPGMMClu_parallel(Meat$x_filtered,
                                                                        Meat$y,
                                                                        g=5,
                                                                        B=B,
                                                                        B.star=B.star))["elapsed"]
# Calculate internal metrics
internal_metrics_UFS <- calculate_internal_metrics(Meat$x,
                                                   out.clu_baseline_UFS$ensemble$label.vec)
# Log the experiment
exp_data <- experiment_logger(
  description = "Clustering with RPGMMClu - UFS:Inf-FS2020",
  dataset = "Meat",
  ensemble_method = "RPGMMClu_parallel",
  ensemble_method_params = list(g = 5, B = B, B.star = B.star),
  execution_time = as.numeric(execution_time),
  labels_clustering = out.clu_baseline_UFS$ensemble$label.vec,
  internal_metrics = internal_metrics_UFS,
  external_metrics = list(ensemble_ari = out.clu_baseline_UFS$ensemble$ari[[1]])
)
save_experiment(exp_data)

experiments_data <- load_experiments()



