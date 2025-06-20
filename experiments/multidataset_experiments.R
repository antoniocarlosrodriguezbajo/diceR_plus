library(diceRplus)
library(mclust)

UFS_Methods <- list(
  "Inf-FS2020" = TRUE
)

algorithms = c("nmf", "hc", "diana", "km", "pam", "ap", "sc",
               "gmm", "block", "som", "cmeans", "hdbscan")


calculate_consensus_labels <- function (E, k, ref_labels) {

  dims <- dim(E)
  E_4d <- array(E, dim = c(dims[1], dims[2], 1, 1))
  dimnames(E_4d) <- list(NULL, NULL, NULL, as.character(k))

  # List of consensus methods and their specific arguments
  consensus_methods <- list(
    kmodes   = function(E) k_modes(E, is.relabelled = FALSE, seed = 1),
    majority = function(E) majority_voting(E, is.relabelled = FALSE),
    lce      = function(E) LCE(E, k, sim.mat = "cts"),
    lca      = function(E) LCA(E, is.relabelled = FALSE, seed = 1),
    cspa     = function(E) CSPA(E_4d, k)
  )

  # Apply each consensus method
  consensus_labels = lapply(consensus_methods, function(f) f(E))

  # Relabel
  aligned_labels = list(
    kmodes   = relabel_class(consensus_labels$kmodes, ref_labels),
    majority = relabel_class(consensus_labels$majority, ref_labels),
    lce      = relabel_class(consensus_labels$lce, ref_labels),
    lca      = relabel_class(consensus_labels$lca, ref_labels),
    cspa     = relabel_class(consensus_labels$cspa, ref_labels)

  )
  return(aligned_labels)
}

# Calcute internal_metrics based on consensus_evaluate of diceR
calculate_internal_metrics <- function(data,cluster_labels) {
  num_labels <- length(cluster_labels)

  # Create structure for consensus_evaluate
  cc_data = array(cluster_labels, dim = c(num_labels, 1, 1, 1))
  row_names = rownames(data)
  dimnames(cc_data) = list(
    row_names,  # Primer nivel de nombres: nombres de las filas de Meat$x
    "R1",       # Repetition
    "RPGMMClu", # Clustering algorithm
    "5"         # Number of clusters
  )
  # Evaluation
  result_evaluation = consensus_evaluate(data, cc_data)
  metrics_df = result_evaluation$ii[[1]]
  metrics_list = lapply(metrics_df[1, -1], function(x) unname(x))
  return(metrics_list)
}

calculate_external_metrics <- function (consensus_labels, ref_labels) {
  # Evaluate results
  eval_results = lapply(consensus_labels, function(labels) list(
    confmat = ev_confmat(labels, ref_labels),
    nmi     = ev_nmi(labels, ref_labels),
    ari     = adjustedRandIndex(labels, ref_labels)
  ))

  return(eval_results)
}

run_RPClu_experiments <- function(data_all, top_features, ref_labels, algorithms, B, k) {
  data <- data_all[,top_features]
  experiment_ids <- c()
  for (alg in algorithms) {
    print(alg)
    tryCatch({
      execution_time = system.time(RPClu_results <- RPClu_parallel(data,
                                                                    clust_fun = alg,
                                                                    g = k,
                                                                    B = B,
                                                                    verb = TRUE))["elapsed"]

      E = as.matrix(RPClu_results$clusterings)

      if (B > 1) {
        consensus_labels = calculate_consensus_labels(E, k, ref_labels)
      } else {
        consensus_labels = list(individual = as.vector(E))
      }

      external_metrics = calculate_external_metrics(consensus_labels, ref_labels)
      for (con_method in names(external_metrics)) {
        ari_valor = round(external_metrics[[con_method]]$ari, 4)
        nmi_valor = round(external_metrics[[con_method]]$nmi, 4)
        acc_valor = external_metrics[[con_method]]$confmat$.estimate[external_metrics[[con_method]]$confmat$.metric == "accuracy"]
        acc_valor = round(acc_valor, 4)
        cat(sprintf("%s: ARI:%.4f NMI:%.4f ACC:%.4f \n",
                    con_method, ari_valor, nmi_valor, acc_valor))

        # Log the experiment
        exp_data = experiment_logger(
          description = paste0("Ensemble Clustering with Inf-FS2020 + RP +",
                               alg,"+", con_method),
          dataset = "Meat",
          clustering_method = alg,
          clustering_method_params = list(k = k, B = B),
          UFS_method = "Inf-FS2020",
          num_features = ncol(data),
          features = top_features,
          ensemble_method = con_method,
          ensemble_method_params = list(k = k, B = B),
          execution_time = as.numeric(execution_time),
          labels_clustering = consensus_labels[[con_method]],
          # Calculate internal metrics
          internal_metrics = calculate_internal_metrics(data,consensus_labels[[con_method]]),
          # Calculate external metrics
          external_metrics = list(ari = ari_valor, nmi = nmi_valor, acc = acc_valor)
        )
        save_experiment(exp_data)
        experiment_ids <- c(experiment_ids, exp_data$experiment_id)
      }

    },
    error = function(e) {
      message("Error: ", conditionMessage(e))
    })
  }
  return(c(first = experiment_ids[1], last = tail(experiment_ids, 1)))
}

run_experiment_individual_algorithm <- function(data, ref_labels, algorithm, k) {
  clustering_labels = consensus_cluster(data, nk = k,
                                             reps = 1, p.item = 1,
                                             algorithms = algorithm, progress = FALSE)
  clustering_labels = relabel_class(clustering_labels, ref_labels)
  clustering_labels = list(individual = clustering_labels)
  external_metrics = calculate_external_metrics(clustering_labels, ref_labels)
  con_method = "individual"
  ari_valor = round(external_metrics[[con_method]]$ari, 4)
  nmi_valor = round(external_metrics[[con_method]]$nmi, 4)
  acc_valor = external_metrics[[con_method]]$confmat$.estimate[external_metrics[[con_method]]$confmat$.metric == "accuracy"]
  acc_valor = round(acc_valor, 4)
  cat(sprintf("%s: ARI:%.4f NMI:%.4f ACC:%.4f \n",
              con_method, ari_valor, nmi_valor, acc_valor))
}

############################
# Meat
############################
data(Meat)

data <- Meat

mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_Meat_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]

sprintf("The %s method selected %d features.", method, length(top_features))

B <- 10
k <- 5

experiments <- run_RPClu_experiments(data=data$x,
                                     top_features,
                                     ref_labels = as.integer(data$y),
                                     algorithms = algorithms,
                                     B=B,
                                     k=k)
print(experiments)

experiments_data_all <- load_experiments()

filter_meat <- experiments_data_all$experiment_id >= experiments[["first"]] &
  experiments_data_all$experiment_id <= experiments[["last"]]

experiments_meat <- experiments_data_all[filter_meat,]

# Extract ARI values
ari_vals <- sapply(experiments_meat$external_metrics, function(x) x$ari)

# Get indices of the top 10 ARI scores
top10_idx <- order(ari_vals, decreasing = TRUE)[1:10]

# Build the summary table
top10_summary <- data.frame(
  clustering_method = experiments_meat$clustering_method[top10_idx],
  ensemble_method = experiments_meat$ensemble_method[top10_idx],
  ari = round(ari_vals[top10_idx], 4)
)

print(top10_summary)







run_experiment_individual_algorithm(data = data$x,
                                    ref_labels = as.integer(data$y),
                                    algorithm = "gmm",
                                    k=k
                                    )

run_experiment_individual_algorithm(data = data$x[, top_features],
                                    ref_labels = as.integer(data$y),
                                    algorithm = "gmm",
                                    k=k
)

run_RPClu_experiments(data = data$x[, top_features],
                      ref_labels = as.integer(data$y),
                      algorithms = "gmm",
                      B=B,
                      k=k)







############################
# ALLAML  TBC
############################
data(ALLAML)

data <- ALLAML

mean(data$x)
sd(data$x)

############################
# leukemia TBC
############################


############################
# lung_cancer
############################
data(lung_cancer)

data <- lung_cancer

mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_lung_cancer_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]
sprintf("The %s method selected %d features.", method, length(top_features))


B <- 50
k <- 5

run_RPClu_experiments(data=data$x[, top_features],
                      ref_labels = as.integer(data$y),
                      algorithms = algorithms,
                      B=B,
                      k=k)















##################
for (alg in algorithms) {
  print(alg)
  tryCatch({
    execution_time <- system.time(RPClu_results <- RPClu_parallel(data$x[, top_features],
                                                                  clust_fun = alg,
                                                                  g = k,
                                                                  B = B,
                                                                  verb = TRUE))["elapsed"]

    E <- as.matrix(RPClu_results$clusterings)
    ref_labels <- as.integer(data$y)

    if (B > 1) {
      consensus_labels <- calculate_consensus_labels(E, k, ref_labels)
    } else {
      consensus_labels <- list(unitary = as.vector(E))
    }

    external_metrics <- calculate_external_metrics(consensus_labels, ref_labels)
    for (con_method in names(external_metrics)) {
      ari_valor <- round(external_metrics[[con_method]]$ari, 4)
      nmi_valor <- round(external_metrics[[con_method]]$nmi, 4)
      acc_valor <- external_metrics[[con_method]]$confmat$.estimate[external_metrics[[con_method]]$confmat$.metric == "accuracy"]
      acc_valor <- round(acc_valor, 4)
      cat(sprintf("%s: ARI:%.4f NMI:%.4f ACC:%.4f \n",
                  con_method, ari_valor, nmi_valor, acc_valor))
    }

    external_metrics$kmodes$nmi
    external_metrics$kmodes$ari
    external_metrics$kmodes$confmat$.estimate[external_metrics$kmodes$confmat$.metric == "accuracy"]
  },
  error = function(e) {
    message("Error: ", conditionMessage(e))
    NA
  })
}
