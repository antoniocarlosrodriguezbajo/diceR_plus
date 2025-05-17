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
                                               out.clu_p_baseline$ensemble$label.vec)
# Log the experiment
exp_data <- experiment_logger(
  description = "Clustering with RPGMMClu - Experiment in Anderlucci's paper",
  dataset = "Meat",
  ensemble_method = "RPGMMClu",
  ensemble_method_params = list(g = 5, B = 100, B.star = 10),
  execution_time = as.numeric(execution_time),
  labels_clustering = out.clu_p_baseline$ensemble$label.vec,
  internal_metrics = internal_metrics,
  external_metrics = list(ensemble_ari = out.clu_p_baseline$ensemble$ari[[1]])
)
save_experiment(exp_data)

############################################################
# Experiment in Anderlucci's paper (parallel implementation)
############################################################



experiments_data <- load_experiments()




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

UFS_Results <- runUFS(Meat$x, UFS_Methods)
save(UFS_Results, file = "experiments/UFS_Meat.RData")

load("experiments/UFS_Meat.RData")

N <- 1000
UFS_Results$Results$"Inf-FS2020"$Result[[1]][1:N]
UFS_Results$ExecutionTimes$"Inf-FS2020"

UFS_Results$Results$Laplacian$Result[[1]][1:N]
UFS_Results$ExecutionTimes$Laplacian



