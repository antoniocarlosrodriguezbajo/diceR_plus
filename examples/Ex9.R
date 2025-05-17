# Anderlucci's meat with UFS integrated into the package
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

UFS_Methods <- list(
  "InfFS" = FALSE,
  "CFS" = FALSE,
  "Laplacian" = TRUE,
  "DGUFS" = TRUE,
  "UFSOL" = TRUE,
  "SPEC" = TRUE,
  "UDFS" = TRUE,
  "RUFS"  = TRUE
)

data(Meat)

UFS_Results <- runUFS(Meat$x, UFS_Methods)

B <- 100
B.star=10

# Original vs Parallel
# Run RPGMMClu baseline
execution_time <- system.time(out.clu_p_baseline <- RPGMMClu(Meat$x,
                                                             Meat$y,
                                                             g=5,
                                                             B=B,
                                                             B.star=B.star,
                                                             verb=TRUE))["elapsed"]


internat_metrics <- calculate_internal_metrics(Meat$x,
                                               out.clu_p_baseline$ensemble$label.vec)

exp_data <- experiment_logger(
  description = "Baseline clustering with RPGMMClu",
  dataset = "Meat",
  ensemble_method = "RPGMMClu",
  ensemble_method_params = list(g = 5, B = 100, B.star = 10),
  execution_time = as.numeric(execution_time),
  labels_clustering = out.clu_p_baseline$ensemble$label.vec,
  internal_metrics = internat_metrics,
  external_metrics = list(ensemble_ari = out.clu_p_baseline$ensemble$ari[[1]])
)

save_experiment(exp_data)

experiments_data <- load_experiments()

