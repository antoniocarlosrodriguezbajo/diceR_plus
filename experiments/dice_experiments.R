library(diceRplus)
library(mclust)
library(R.matlab)

algorithms = c("nmf", "hc", "diana", "km", "pam", "ap", "sc",
               "gmm", "block", "som", "cmeans", "hdbscan")

UFS_Methods <- list(
  "Inf-FS2020" = TRUE
)


mean_abs_correlation <- function(data) {
  cor_mat <- cor(data, use = "pairwise.complete.obs")
  round(mean(abs(cor_mat[upper.tri(data)])),4)
}

get_best_clustering <- function(data, nk, algorithms,
                                cons.funs = c("kmodes", "majority", "CSPA", "LCE", "LCA"),
                                reps=10) {
  dice.obj <- dice(
    data = data,
    nk = nk,
    algorithms = algorithms,
    cons.funs = cons.funs,
    reps = 10,
    trim = TRUE,
    n = 1
  )

  best_clustering_alg = dice.obj$indices$trim$alg.keep

  cc_data_list <- list()

  for (j in colnames(dice.obj$clusters)) {
    column <- dice.obj$clusters[, j]
    print(paste("Method:", j))
    num_labels = length(column)
    cc_data = array(column, dim = c(num_labels, 1, 1, 1))
    row_names = rownames(column)
    dimnames(cc_data) = list(
      row_names,             # Row names
      "R1",                  # Repetition (just a placeholder)
      j,                     # Consensus clustering algorithm
      as.character(nk)       # Number of clusters
    )
    cc_data_list[[j]] = cc_data
  }

  result = do.call(consensus_evaluate, c(list(data), cc_data_list, list(n = 1, trim = TRUE)))

  best_consensus_alg = result$trim.obj$alg.keep

  return(list(
    best_clustering_alg = best_clustering_alg,
    best_consensus_alg = best_consensus_alg,
    best_clustering = as.vector(result$trim.obj$E.new[[1]])
  ))
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
table(data$y)
mean_abs_correlation(data$x)




result <- get_best_clustering(data = data$x,
                              nk = 5,
                              algorithms = algorithms)

paste(result$best_clustering_alg, result$best_consensus_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))

# "SOM-LCE"
# 0.2440049


# Start MATLAB server
Matlab$startServer()
matlab <- Matlab()
print(matlab)

# Check if MATLAB server is running
isOpen <- open(matlab)
if (!isOpen) {
  print("MATLAB server is not running: waited 30 seconds.")
}
print(matlab)

# Set the dataset variable
setVariable(matlab, X = data$x)

# Evaluate MATLAB command to display variables
evaluate(matlab, "whos")

evaluate(matlab, "rng('default');")
evaluate(matlab, "rng(42);")

for (alpha in seq(0, 1, by = 0.1)) {

  evaluate(matlab, "clear RANKED WEIGHT SUBSET;")


  # Run the feature selection method in MATLAB
  cmd <- sprintf('[RANKED, WEIGHT, SUBSET] = InfFS_U(X, %.1f);', alpha)
  evaluate(matlab, cmd)

  ranked  <- getVariable(matlab, "RANKED")
  weight  <- getVariable(matlab, "WEIGHT")
  subset  <- getVariable(matlab, "SUBSET")

  top_features <- as.vector(subset$SUBSET)
  num_features = length(top_features)

  result <- get_best_clustering(data = data$x[,top_features],
                                nk = 5,
                                algorithms = "som",
                                cons.funs = "LCE")
  print(paste(
    alpha,
    num_features,
    result$best_clustering_alg, result$best_consensus_alg,
    sep = "-"))

  print(adjustedRandIndex(result$best_clustering, as.numeric(data$y)))

}

# Close MATLAB connection
close(matlab)


# [1] "Method: LCE"
# [1] "0.1-91-SOM-LCE"
# [1] 0.3821559

########################################################
########################################################
# Lymphoma
########################################################
########################################################

data(lymphoma)

data <- lymphoma

mean(data$x)
sd(data$x)
table(data$y)
mean_abs_correlation(data$x)




result <- get_best_clustering(data = data$x,
                              nk = 2,
                              algorithms = algorithms)

paste(result$best_clustering_alg, result$best_consensus_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))

# > paste(result$best_clustering_alg, result$best_consensus_alg,  sep = "-")
# [1] "SC-CSPA"
# > adjustedRandIndex(result$best_clustering, as.numeric(data$y))
# [1] 0.8306733


# Start MATLAB server
Matlab$startServer()
matlab <- Matlab()
print(matlab)

# Check if MATLAB server is running
isOpen <- open(matlab)
if (!isOpen) {
  print("MATLAB server is not running: waited 30 seconds.")
}
print(matlab)

# Set the dataset variable
setVariable(matlab, X = data$x)

# Evaluate MATLAB command to display variables
evaluate(matlab, "whos")

evaluate(matlab, "rng('default');")
evaluate(matlab, "rng(42);")

for (alpha in seq(0, 1, by = 0.1)) {

  evaluate(matlab, "clear RANKED WEIGHT SUBSET;")


  # Run the feature selection method in MATLAB
  cmd <- sprintf('[RANKED, WEIGHT, SUBSET] = InfFS_U(X, %.1f);', alpha)
  evaluate(matlab, cmd)

  ranked  <- getVariable(matlab, "RANKED")
  weight  <- getVariable(matlab, "WEIGHT")
  subset  <- getVariable(matlab, "SUBSET")

  top_features <- as.vector(subset$SUBSET)
  num_features = length(top_features)

  result <- get_best_clustering(data = data$x[,top_features],
                                nk = 2,
                                algorithms = "sc",
                                cons.funs = "CSPA")
  print(paste(
    alpha,
    num_features,
    result$best_clustering_alg, result$best_consensus_alg,
    sep = "-"))

  print(adjustedRandIndex(result$best_clustering, as.numeric(data$y)))

}

# Close MATLAB connection
close(matlab)

