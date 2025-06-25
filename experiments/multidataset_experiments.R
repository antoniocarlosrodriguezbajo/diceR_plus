library(diceRplus)
library(mclust)
library(scales)

UFS_Methods <- list(
  "Inf-FS2020" = TRUE
)

algorithms = c("nmf", "hc", "diana", "km", "pam", "ap", "sc",
               "gmm", "block", "som", "cmeans", "hdbscan")

algorithms = c("gmm")


calculate_consensus_labels <- function(E, k, ref_labels, method = NULL) {
  # Reshape the input matrix E into a 4D array for methods that require it (e.g., CSPA)
  dims <- dim(E)
  E_4d <- array(E, dim = c(dims[1], dims[2], 1, 1))
  dimnames(E_4d) <- list(NULL, NULL, NULL, as.character(k))

  # Define the list of available consensus methods and their parameters
  consensus_methods <- list(
    #kmodes   = function(E) k_modes(E, is.relabelled = FALSE, seed = 1),
    #majority = function(E) majority_voting(E, is.relabelled = FALSE),
    LCE      = function(E) LCE(E, k, sim.mat = "cts")
    #LCA      = function(E) LCA(E, is.relabelled = FALSE, seed = 1),
    #CSPA     = function(E) CSPA(E_4d, k)
  )

  # If a specific method is provided, apply only that one
  if (!is.null(method)) {
    if (!method %in% names(consensus_methods)) {
      stop("Unrecognized method: ", method)
    }
    consensus_labels <- list()
    consensus_labels[[method]] <- consensus_methods[[method]](E)
  } else {
    # Otherwise, apply all available methods
    consensus_labels <- lapply(consensus_methods, function(f) f(E))
  }

  # Relabel the consensus results to match the reference labels
  aligned_labels <- lapply(names(consensus_labels), function(m) {
    relabel_class(consensus_labels[[m]], ref_labels)
  })
  names(aligned_labels) <- names(consensus_labels)

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
    "R1",       # Repetition (just a placeholder)
    "RPGMMClu", # Clustering algorithm (just a placeholder)
    "5"         # Number of clusters (just a placeholder)
  )
  # Evaluation
  result_evaluation = consensus_evaluate(data, cc_data)
  metrics_df = result_evaluation$ii[[1]]
  metrics_list = lapply(metrics_df[1, -1], function(x) unname(x))
  return(metrics_list)
}

calculate_external_metrics <- function (consensus_labels, ref_labels, method = NULL) {
  # Evaluate results
  eval_results = lapply(consensus_labels, function(labels) list(
    confmat = ev_confmat(labels, ref_labels),
    nmi     = ev_nmi(labels, ref_labels),
    ari     = adjustedRandIndex(labels, ref_labels)
  ))

  return(eval_results)
}

run_RPClu_experiments <- function(data_all,
                                  top_features,
                                  ref_labels,
                                  algorithms,
                                  B, k,
                                  dataset = "Meat",
                                  UFS_method = "Inf-FS2020",
                                  consensus_method = NULL,
                                  rp = TRUE,
                                  method_rp = "gaussian",
                                  seed = 101) {
  if (is.null(top_features)) {
    data = data_all
  } else {
    data = data_all[,top_features]
  }
  experiment_ids = c()
  for (alg in algorithms) {
    print(alg)
    tryCatch({
      execution_time <- system.time({
        RPClu_results <- RPClu_parallel(data,
                                        clust_fun = alg,
                                        g = k,
                                        B = B,
                                        rp = rp,
                                        method_rp = method_rp,
                                        seed = seed)
      })["elapsed"]


      E = as.matrix(RPClu_results$clusterings)

      if (B > 1) {
        consensus_labels = calculate_consensus_labels(E, k, ref_labels, consensus_method)
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
        if (rp) {
          description = paste0("Ensemble Clustering + ", UFS_method, " + RP + ",
                               alg, " + " , con_method, " + " ,
                               " seed ", as.character(seed))
        } else {
          description = paste0("Ensemble Clustering + ", UFS_method, " + NO RP + ",
                               alg, " + " , con_method, " + " ,
                               " seed ", as.character(seed))
        }
        exp_data = experiment_logger(
          description = description,
          dataset = dataset,
          clustering_method = alg,
          clustering_method_params = list(k = k, B = B),
          UFS_method = UFS_method,
          num_features = ncol(data),
          features = top_features,
          consensus_method = con_method,
          consensus_method_params = list(k = k, B = B),
          execution_time = as.numeric(execution_time),
          labels_clustering = consensus_labels[[con_method]],
          # Calculate internal metrics
          internal_metrics = calculate_internal_metrics(data,consensus_labels[[con_method]]),
          # Calculate external metrics
          external_metrics = list(ari = ari_valor, nmi = nmi_valor, acc = acc_valor)
        )
        save_experiment(exp_data)
        experiment_ids = c(experiment_ids, exp_data$experiment_id)
      }

    },
    error = function(e) {
      message("Error: ", conditionMessage(e))
    })
  }
  return(list(first = experiment_ids[1], last = tail(experiment_ids, 1)))
}


run_clustering_experiments <- function(data_all,
                                       top_features,
                                       ref_labels,
                                       algorithms,
                                       k,
                                       dataset = "Meat",
                                       UFS_method = "Inf-FS2020",
                                       consensus_method = NULL,
                                       seed = 101) {
  if (is.null(top_features)) {
    data = data_all
  } else {
    data = data_all[,top_features]
  }
  set.seed(seed)
  execution_time = system.time({
    results_dice = dice(data = data,
                         nk = k,
                         algorithms = algorithms,
                         cons.funs = consensus_method,
                         evaluate = TRUE,
                         ref.cl = ref_labels,
                         seed = seed,
                         seed.data = seed,
                         progress = FALSE,
                         verbose = FALSE)
  })["elapsed"]

  consensus_labels = as.vector(results_dice$clusters[ , 2])
  # Calculate internal metrics
  internal_metrics = calculate_internal_metrics(data,consensus_labels)
  # Calculate external metrics
  ei_k = results_dice$indices$ei[[as.character(k)]]
  ei_filtered = ei_k[ei_k$Algorithms == consensus_method, -1]
  acc_value = round(as.numeric(ei_filtered[["accuracy"]]),4)
  nmi_value = round(as.numeric(ei_filtered[["NMI"]]),4)
  ari_value = round(adjustedRandIndex(consensus_labels, ref_labels),4)

  description = paste0("Ensemble Clustering + ", UFS_method, " + NO RP + ",
                       algorithms, " + " , consensus_method, " + " ,
                       " seed ", as.character(seed))
  exp_data = experiment_logger(
    description = description,
    dataset = dataset,
    clustering_method = algorithms,
    clustering_method_params = list(k = k),
    UFS_method = UFS_method,
    num_features = ncol(data),
    features = top_features,
    consensus_method = consensus_method,
    consensus_method_params = list(k = k),
    execution_time = as.numeric(execution_time),
    labels_clustering = consensus_labels,
    internal_metrics = internal_metrics,
    external_metrics = list(ari = ari_value, nmi = nmi_value, acc = acc_value)
  )
  save_experiment(exp_data)
  return(exp_data$experiment_id)
}


summarize_top_metrics <- function(experiments, top_n = 10) {
  # Load all saved experiments
  experiments_data_all <- load_experiments()

  # Filter to include only the experiments in the relevant ID range
  experiment_filter <- experiments_data_all$experiment_id >= experiments$first &
    experiments_data_all$experiment_id <= experiments$last
  experiments_dataset <- experiments_data_all[experiment_filter, ]

  # Extract external evaluation metrics
  ari_vals <- sapply(experiments_dataset$external_metrics, function(x) x$ari)
  nmi_vals <- sapply(experiments_dataset$external_metrics, function(x) x$nmi)
  acc_vals <- sapply(experiments_dataset$external_metrics, function(x) x$acc)

  # Get indices of top-N values for each metric
  top_ari_idx <- order(ari_vals, decreasing = TRUE)[1:top_n]
  top_nmi_idx <- order(nmi_vals, decreasing = TRUE)[1:top_n]
  top_acc_idx <- order(acc_vals, decreasing = TRUE)[1:top_n]

  # Build summary tables for each metric
  top_ari <- data.frame(
    metric = "ARI",
    clustering_method = experiments_dataset$clustering_method[top_ari_idx],
    consensus_method = experiments_dataset$consensus_method[top_ari_idx],
    value = round(ari_vals[top_ari_idx], 4)
  )

  top_nmi <- data.frame(
    metric = "NMI",
    clustering_method = experiments_dataset$clustering_method[top_nmi_idx],
    consensus_method = experiments_dataset$consensus_method[top_nmi_idx],
    value = round(nmi_vals[top_nmi_idx], 4)
  )

  top_acc <- data.frame(
    metric = "ACC",
    clustering_method = experiments_dataset$clustering_method[top_acc_idx],
    consensus_method = experiments_dataset$consensus_method[top_acc_idx],
    value = round(acc_vals[top_acc_idx], 4)
  )

  # Combine all into one summary table
  top_summary <- rbind(top_ari, top_nmi, top_acc)
  return(top_summary)
}

evaluate_experiments <- function(first_id, n_reps) {
  last_id <- first_id + n_reps - 1

  all_experiments <- load_experiments()

  filter_experiments <- all_experiments$experiment_id >= first_id &
    all_experiments$experiment_id <= last_id

  experiments_subset <- all_experiments[filter_experiments, ]

  extract_metric <- function(metric_name) {
    vals <- sapply(experiments_subset$external_metrics, function(x) x[[metric_name]])
    list(mean = mean(vals, na.rm = TRUE), sd = sd(vals, na.rm = TRUE))
  }

  ari <- extract_metric("ari")
  nmi <- extract_metric("nmi")
  acc <- extract_metric("acc")

  cat("Mean ARI:", round(ari$mean, 4), "| SD ARI:", round(ari$sd, 4), "\n")
  cat("Mean NMI:", round(nmi$mean, 4), "| SD NMI:", round(nmi$sd, 4), "\n")
  cat("Mean ACC:", round(acc$mean, 4), "| SD ACC:", round(acc$sd, 4), "\n")
}

run_pipeline <- function(dataset, dataset_name, k, top_features, top_metric="ARI") {
  B <- ceiling(length(top_features) / 100) * 10
  B <- min(B, 100)

  print(sprintf("%d rows, %d features, %d top features, %d random projections, %d clusters",
                nrow(dataset$x), ncol(dataset$x), length(top_features), B, k))

  print("Selecting clustering and consensus methods")
  seed <- 101
  experiments <- run_RPClu_experiments(data_all=dataset$x,
                                       top_features = top_features,
                                       ref_labels = as.integer(dataset$y),
                                       algorithms = algorithms,
                                       B=B, k=k,
                                       dataset = dataset_name,
                                       UFS_method = "Inf-FS2020",
                                       rp = TRUE,
                                       seed = seed)
  print(experiments)

  top_summary <- summarize_top_metrics(experiments, top_n = 10)
  top_metric <- subset(top_summary, metric == top_metric)
  top_metric <- top_metric[order(-top_metric$value, top_metric$clustering_method), ]
  print(top_metric)

  best_clustering <- top_metric$clustering_method[1]
  best_consensus <- top_metric$consensus_method[1]

  print(sprintf("Best Clustering: %s, Best Consensus: %s",
                best_clustering, best_consensus))


  ## Repetitions of best clustering and best consensus (RP + UFS)

  n_reps <- 20
  for (i in 1:n_reps) {
    seed <- 100 + i
    experiments_rep <- run_RPClu_experiments(data_all = dataset$x,
                                             top_features = top_features,
                                             ref_labels = as.integer(dataset$y),
                                             algorithms = best_clustering,
                                             consensus_method = best_consensus,
                                             B = B, k = k,
                                             dataset = dataset_name,
                                             UFS_method = "Inf-FS2020",
                                             rp = TRUE,
                                             seed = seed)
    print(experiments_rep$first)
  }
  print("Repetitions of best clustering and best consensus (RP + UFS)")
  evaluate_experiments(experiments$last + 1,
                       experiments$last + 1 + n_reps - 1)


  experiments$last <- experiments$last + 1 + n_reps - 1


  ## Repetitions of best clustering and best consensus (RP + NO UFS)

  n_reps <- 20
  for (i in 1:n_reps) {
    seed <- 100 + i
    experiments_rep <- run_RPClu_experiments(data_all = dataset$x,
                                             top_features = NULL,
                                             ref_labels = as.integer(dataset$y),
                                             algorithms = best_clustering,
                                             consensus_method = best_consensus,
                                             B = B, k = k,
                                             dataset = dataset_name,
                                             UFS_method = "No UFS",
                                             rp = TRUE,
                                             seed = seed)
    print(experiments_rep$first)
  }

  print("Repetitions of best clustering and best consensus (RP + NO UFS)")
  evaluate_experiments(experiments$last + 1,
                       experiments$last + 1 + n_reps - 1)

  experiments$last <- experiments$last + 1 + n_reps - 1

  ## Repetitions of best clustering and best consensus (NO RP + UFS)

  n_reps <- 20
  for (i in 1:n_reps) {
    seed <- 100 + i
    experiment_id <- run_clustering_experiments(data_all = dataset$x,
                                                top_features = top_features,
                                                ref_labels = as.integer(dataset$y),
                                                algorithms = best_clustering,
                                                consensus_method = best_consensus,
                                                k = k,
                                                dataset = dataset_name,
                                                UFS_method = "Inf-FS2020",
                                                seed = seed)
    print(experiment_id)
  }

  print("Repetitions of best clustering and best consensus (NO RP + UFS)")
  evaluate_experiments(experiments$last + 1,
                       experiments$last + 1 + n_reps - 1)

  experiments$last <- experiments$last + 1 + n_reps - 1

  ## Repetitions of best clustering and best consensus (NO RP + NO UFS)

  n_reps <- 20
  for (i in 1:n_reps) {
    seed <- 100 + i
    experiment_id <- run_clustering_experiments(data_all = dataset$x,
                                                top_features = NULL,
                                                ref_labels = as.integer(dataset$y),
                                                algorithms = best_clustering,
                                                consensus_method = best_consensus,
                                                k = k,
                                                dataset = dataset_name,
                                                UFS_method = "No UFS",
                                                seed = seed)
    print(experiment_id)
  }

  print("Repetitions of best clustering and best consensus (NO RP + NO UFS)")
  evaluate_experiments(experiments$last + 1,
                       experiments$last + 1 + n_reps - 1)

  experiments$last <- experiments$last + 1 + n_reps - 1


}


########################################################
########################################################
# Meat
########################################################
########################################################

data(Meat)

data <- Meat

mean(data$x)
sd(data$x)

# UFS
UFS_results_file <- "experiments/UFS_Meat_Inf-FS2020.RData"
UFS_results <- runUFS(data$x, UFS_Methods)
save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]

run_pipeline(dataset = data,
             dataset_name = "Meat",
             k = 5,
             top_features = top_features,
             top_metric = "ARI")



########################################################
########################################################
# lung_cancer
########################################################
########################################################

data(lung_cancer)

data <- lung_cancer

mean(data$x)
sd(data$x)

# UFS
UFS_results_file <- "experiments/UFS_lung_cancer_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]

run_pipeline(dataset = data,
             dataset_name = "Lung Cancer",
             k = 5,
             top_features = top_features,
             top_metric = "ARI")


########################################################
########################################################
# lymphoma
########################################################
########################################################

data(lymphoma)

data <- lymphoma

data$y <- data$y+1

mean(data$x)
sd(data$x)

# UFS
UFS_results_file <- "experiments/UFS_lymphoma_Inf-FS2020.RData"
UFS_results <- runUFS(data$x, UFS_Methods)
save(UFS_results, file = UFS_results_file)

load(UFS_results_file)
method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]

run_pipeline(dataset = data,
             dataset_name = "lymphoma",
             k = 3,
             top_features = top_features,
             top_metric = "ARI")



########################################################
########################################################
# prostate_ge
########################################################
########################################################

data(prostate_ge)

data <- prostate_ge

mean(data$x)
sd(data$x)

# UFS
UFS_results_file <- "experiments/UFS_prostate_ge_Inf-FS2020.RData"
UFS_results <- runUFS(data$x, UFS_Methods)
save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]

run_pipeline(dataset = data,
             dataset_name = "prostate_ge",
             k = 2,
             top_features = top_features,
             top_metric = "ACC")


########################################################
########################################################
# warpPIE10P
########################################################
########################################################

data(warpPIE10P)

data <- warpPIE10P

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

data$x <- apply(data$x, 2, rescale)

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

# UFS
UFS_results_file <- "experiments/UFS_warpPIE10P_Inf-FS2020.RData"
UFS_results <- runUFS(data$x, UFS_Methods)
save(UFS_results, file = UFS_results_file)

load(UFS_results_file)
method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]

run_pipeline(dataset = data,
             dataset_name = "warpPIE10P",
             k = 10,
             top_features = top_features,
             top_metric = "ACC")




########################################################
########################################################
# ALLAML
########################################################
########################################################

data(ALLAML)

data <- ALLAML

mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_ALLAML_Inf-FS2020.RData"
UFS_results <- runUFS(data$x, UFS_Methods)
save(UFS_results, file = UFS_results_file)

load(UFS_results_file)
method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]

run_pipeline(dataset = data,
             dataset_name = "ALLAML",
             k = 2,
             top_features = top_features,
             top_metric = "ACC")
