% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experiment_logger.R
\name{save_experiment}
\alias{save_experiment}
\title{Save experiment results to an RDS file}
\usage{
save_experiment(experiment_data, file = "experiments/experiment_results.rds")
}
\arguments{
\item{experiment_data}{A dataframe containing the experiment results.}

\item{file}{A character string specifying the file path for storing the experiments.
Default is "experiments/experiment_results.rds".}
}
\value{
NULL (saves the updated experiment log to disk).
}
\description{
This function stores experiment data in an RDS file. If the file already exists,
it appends the new experiment data to the existing records.
}
