#' Generate a unique experiment ID
#'
#' This function retrieves the last experiment ID from an RDS file and increments it.
#'
#' @param file A character string specifying the file path for retrieving experiments.
#' Default is "experiments/experiment_results.rds".
#'
#' @return Integer representing the next available experiment ID.
#' @export
generate_experiment_id <- function(file = "experiments/experiment_results.rds") {
  if (!file.exists(file)) {
    return(1)  # First experiment
  } else {
    experiments <- readRDS(file)
    return(max(experiments$experiment_id, na.rm = TRUE) + 1)  # Increment ID
  }
}

#' Log an experiment result
#'
#' @param description Character, description of experiment.
#' @param dataset Character, name of the dataset.
#' @param ensemble_method Character, ensemble clustering method.
#' @param ensemble_method_params List, parameters of the ensemble clustering method.
#' @param UFS_method Character, unsupervised feature selection method.
#' @param UFS_method_params List, parameters of the unsupervised feature selection method.
#' @param num_features Integer, number of features selected.
#' @param features List, features selected.
#' @param dim_reduction_method Character, dimensionality reduction method.
#' @param dim_reduction_method_params List, parameters of the dimensionality reduction method.
#' @param execution_time Numeric, execution time in seconds.
#' @param labels_clustering ensemble clustering labels.
#' @param internal_metrics List, internal evaluation metrics.
#' @param external_metrics List, external evaluation metrics.
#'
#' @return Dataframe containing the experiment results.
#' @export
experiment_logger <- function(description,
                                   dataset,
                                   clustering_method,
                                   clustering_method_params = list(),
                                   ensemble_method = NA,
                                   ensemble_method_params = list(),
                                   UFS_method = "N/A",
                                   UFS_method_params = list(),
                                   num_features = NA,
                                   features = list(),
                                   dim_reduction_method = "N/A",
                                   dim_reduction_method_params = list(),
                                   execution_time = NA,
                                   labels_clustering = list(),
                                   internal_metrics = list(),
                                   external_metrics = list()) {
  data.frame(
    experiment_id = generate_experiment_id(),
    description = description,
    dataset = dataset,
    clustering_method = clustering_method,
    clustering_method_params = I(list(clustering_method_params)),
    ensemble_method = ensemble_method,
    ensemble_method_params = I(list(ensemble_method_params)),
    UFS_method = UFS_method,
    UFS_method_params = I(list(UFS_method_params)),
    num_features = num_features,
    features = I(list(features)),
    dim_reduction_method = dim_reduction_method,
    dim_reduction_method_params = I(list(dim_reduction_method_params)),
    execution_time = execution_time,
    timestamp = Sys.time(),
    labels_clustering = I(list(labels_clustering)),
    internal_metrics = I(list(internal_metrics)),
    external_metrics = I(list(external_metrics)),
    stringsAsFactors = FALSE
  )
}

#' Save experiment results to an RDS file
#'
#' This function stores experiment data in an RDS file. If the file already exists,
#' it appends the new experiment data to the existing records.
#'
#' @param experiment_data A dataframe containing the experiment results.
#' @param file A character string specifying the file path for storing the experiments.
#' Default is "experiments/experiment_results.rds".
#'
#' @return NULL (saves the updated experiment log to disk).
#' @export
save_experiment <- function(experiment_data, file = "experiments/experiment_results.rds") {
  if (!file.exists(file)) {
    experiments <- experiment_data  # First experiment
  } else {
    experiments <- readRDS(file)  # Load existing experiments
    experiments <- rbind(experiments, experiment_data)  # Append new row
  }

  saveRDS(experiments, file)  # Save updated log
}

#' Load experiment results from an RDS file
#'
#' @param file A character string specifying the file path for retrieving experiments.
#' Default is "experiments/experiment_results.rds".
#'
#' @return Dataframe containing all stored experiments.
#' @export
load_experiments <- function(file = "experiments/experiment_results.rds") {
  if (file.exists(file)) {
    return(readRDS(file))
  } else {
    message("No experiment log found.")
    return(NULL)
  }
}

