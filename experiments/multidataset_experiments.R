library(diceRplus)
UFS_Methods <- list(
  "Inf-FS2020" = TRUE
)

algorithms = c("nmf", "hc", "diana", "km", "pam", "ap", "sc",
               "gmm", "block", "som", "cmeans", "hdbscan")

calculate_external_metrics <- function (E, k, ref_labels) {
  # List of consensus methods and their specific arguments
  consensus_methods <- list(
    kmodes   = function(E) k_modes(E, is.relabelled = FALSE, seed = 1),
    majority = function(E) majority_voting(E, is.relabelled = FALSE),
    lce      = function(E) LCE(E, k, sim.mat = "cts"),
    lca      = function(E) LCA(E, is.relabelled = FALSE, seed = 1)
  )

  # Apply each consensus method
  consensus_labels <- lapply(consensus_methods, function(f) f(E))

  # Relabel only the methods that need it (adjust as needed)
  aligned_labels <- list(
    kmodes   = relabel_class(consensus_labels$kmodes, ref_labels),
    majority = relabel_class(consensus_labels$majority, ref_labels),
    lce      = relabel_class(consensus_labels$lce, ref_labels),
    lca      = relabel_class(consensus_labels$lca, ref_labels)
  )

  # Evaluate results
  eval_results <- lapply(aligned_labels, function(labels) list(
    confmat = ev_confmat(labels, ref_labels),
    nmi     = ev_nmi(labels, ref_labels),
    ari     = adjustedRandIndex(labels, ref_labels)
  ))

  return(eval_results)
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

B <- 100
execution_time <- system.time(RPClu_results <- RPClu_parallel(data$x[, top_features],
                                                                 clust_fun = "gmm",
                                                                 g=5,
                                                                 B=B,
                                                                 verb=TRUE))["elapsed"]

E <- as.matrix(RPClu_results$clusterings)

k=5
ref_labels <- as.integer(data$y)

external_metrics <- calculate_external_metrics(E,k,ref_labels)

external_metrics

# > external_metrics
# $kmodes
# $kmodes$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.823
# 2 kap                  multiclass     0.775
# 3 sens                 macro          0.843
# 4 spec                 macro          0.954
# 5 ppv                  macro          0.847
# 6 npv                  macro          0.954
# 7 mcc                  multiclass     0.776
# 8 j_index              macro          0.797
# 9 bal_accuracy         macro          0.898
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.847
# 12 recall               macro          0.843
# 13 f_meas               macro          0.843
#
# $kmodes$nmi
# [1] 0.7794098
#
# $kmodes$ari
# [1] 0.666964
#
#
# $majority
# $majority$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.814
# 2 kap                  multiclass     0.764
# 3 sens                 macro          0.836
# 4 spec                 macro          0.951
# 5 ppv                  macro          0.872
# 6 npv                  macro          0.958
# 7 mcc                  multiclass     0.786
# 8 j_index              macro          0.787
# 9 bal_accuracy         macro          0.894
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.872
# 12 recall               macro          0.836
# 13 f_meas               macro          0.821
#
# $majority$nmi
# [1] 0.8084113
#
# $majority$ari
# [1] 0.6827528
#
#
# $lce
# $lce$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.745
# 2 kap                  multiclass     0.676
# 3 sens                 macro          0.778
# 4 spec                 macro          0.933
# 5 ppv                  macro          0.684
# 6 npv                  macro          0.948
# 7 mcc                  multiclass     0.730
# 8 j_index              macro          0.711
# 9 bal_accuracy         macro          0.855
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.684
# 12 recall               macro          0.778
# 13 f_meas               macro          0.713
#
# $lce$nmi
# [1] 0.8507881
#
# $lce$ari
# [1] 0.6882715
#
#
# $lca
# $lca$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.805
# 2 kap                  multiclass     0.753
# 3 sens                 macro          0.829
# 4 spec                 macro          0.949
# 5 ppv                  macro          0.837
# 6 npv                  macro          0.951
# 7 mcc                  multiclass     0.759
# 8 j_index              macro          0.778
# 9 bal_accuracy         macro          0.889
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.837
# 12 recall               macro          0.829
# 13 f_meas               macro          0.825
#
# $lca$nmi
# [1] 0.7791981
#
# $lca$ari
# [1] 0.6592031
#
# "The Inf-FS2020 method selected 91 features."

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

g <- 5
B <- 100
execution_time <- system.time(RPClu_results <- RPClu_parallel(data$x[, top_features],
                                                              clust_fun = "gmm",
                                                              g=g,
                                                              B=B,
                                                              verb=TRUE))["elapsed"]

E <- as.matrix(RPClu_results$clusterings)

external_metrics <- calculate_external_metrics(E,g, as.integer(data$y))

external_metrics

# > external_metrics
# $kmodes
# $kmodes$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.690
# 2 kap                  multiclass     0.497
# 3 sens                 macro          0.634
# 4 spec                 macro          0.912
# 5 ppv                  macro          0.646
# 6 npv                  macro          0.888
# 7 mcc                  multiclass     0.527
# 8 j_index              macro          0.546
# 9 bal_accuracy         macro          0.773
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.646
# 12 recall               macro          0.634
# 13 f_meas               macro          0.627
#
# $kmodes$nmi
# [1] 0.4651067
#
# $kmodes$ari
# [1] 0.3685128
#
#
# $majority
# $majority$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.704
# 2 kap                  multiclass     0.499
# 3 sens                 macro          0.621
# 4 spec                 macro          0.911
# 5 ppv                  macro          0.696
# 6 npv                  macro          0.890
# 7 mcc                  multiclass     0.524
# 8 j_index              macro          0.532
# 9 bal_accuracy         macro          0.766
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.696
# 12 recall               macro          0.621
# 13 f_meas               macro          0.634
#
# $majority$nmi
# [1] 0.4598495
#
# $majority$ari
# [1] 0.4061217
#
#
# $lce
# $lce$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.700
# 2 kap                  multiclass     0.483
# 3 sens                 macro          0.577
# 4 spec                 macro          0.912
# 5 ppv                  macro          0.588
# 6 npv                  macro          0.890
# 7 mcc                  multiclass     0.504
# 8 j_index              macro          0.489
# 9 bal_accuracy         macro          0.744
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.588
# 12 recall               macro          0.577
# 13 f_meas               macro          0.560
#
# $lce$nmi
# [1] 0.4692405
#
# $lce$ari
# [1] 0.4093116
#
#
# $lca
# $lca$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.611
# 2 kap                  multiclass     0.429
# 3 sens                 macro          0.657
# 4 spec                 macro          0.902
# 5 ppv                  macro          0.579
# 6 npv                  macro          0.877
# 7 mcc                  multiclass     0.480
# 8 j_index              macro          0.559
# 9 bal_accuracy         macro          0.779
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.579
# 12 recall               macro          0.657
# 13 f_meas               macro          0.579
#
# $lca$nmi
# [1] 0.4035596
#
# $lca$ari
# [1] 0.2510463
#
# [1] "The Inf-FS2020 method selected 583 features."

############################
# lymphoma TBC
############################
data(lymphoma)

data <- lymphoma

mean(data$x)
sd(data$x)

############################
# prostate_ge
############################

data(prostate_ge)

data <- prostate_ge

mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_prostate_ge_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]
sprintf("The %s method selected %d features.", method, length(top_features))

g <- 2
B <- 100
execution_time <- system.time(RPClu_results <- RPClu_parallel(data$x[, top_features],
                                                              clust_fun = "km",
                                                              g=g,
                                                              B=B,
                                                              verb=TRUE))["elapsed"]

E <- as.matrix(RPClu_results$clusterings)

external_metrics <- calculate_external_metrics(E,g, as.integer(data$y))

external_metrics

# > external_metrics
# $kmodes
# $kmodes$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             binary         0.716
# 2 kap                  binary         0.429
# 3 sens                 binary         0.62
# 4 spec                 binary         0.808
# 5 ppv                  binary         0.756
# 6 npv                  binary         0.689
# 7 mcc                  binary         0.436
# 8 j_index              binary         0.428
# 9 bal_accuracy         binary         0.714
# 10 detection_prevalence binary         0.402
# 11 precision            binary         0.756
# 12 recall               binary         0.62
# 13 f_meas               binary         0.681
#
# $kmodes$nmi
# [1] 0.1444463
#
# $kmodes$ari
# [1] 0.1782498
#
#
# $majority
# $majority$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             binary         0.569
# 2 kap                  binary         0.127
# 3 sens                 binary         0.26
# 4 spec                 binary         0.865
# 5 ppv                  binary         0.65
# 6 npv                  binary         0.549
# 7 mcc                  binary         0.158
# 8 j_index              binary         0.125
# 9 bal_accuracy         binary         0.563
# 10 detection_prevalence binary         0.196
# 11 precision            binary         0.65
# 12 recall               binary         0.26
# 13 f_meas               binary         0.371
#
# $majority$nmi
# [1] 0.02151663
#
# $majority$ari
# [1] 0.01253643
#
#
# $lce
# $lce$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             binary        0.569
# 2 kap                  binary        0.123
# 3 sens                 binary        0.14
# 4 spec                 binary        0.981
# 5 ppv                  binary        0.875
# 6 npv                  binary        0.543
# 7 mcc                  binary        0.225
# 8 j_index              binary        0.121
# 9 bal_accuracy         binary        0.560
# 10 detection_prevalence binary        0.0784
# 11 precision            binary        0.875
# 12 recall               binary        0.14
# 13 f_meas               binary        0.241
#
# $lce$nmi
# [1] 0.06406633
#
# $lce$ari
# [1] 0.01575351
#
#
# $lca
# $lca$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             binary         0.755
# 2 kap                  binary         0.509
# 3 sens                 binary         0.72
# 4 spec                 binary         0.788
# 5 ppv                  binary         0.766
# 6 npv                  binary         0.745
# 7 mcc                  binary         0.510
# 8 j_index              binary         0.508
# 9 bal_accuracy         binary         0.754
# 10 detection_prevalence binary         0.461
# 11 precision            binary         0.766
# 12 recall               binary         0.72
# 13 f_meas               binary         0.742
#
# $lca$nmi
# [1] 0.1971804
#
# $lca$ari
# [1] 0.2525461
#
# ] "The Inf-FS2020 method selected 949 features."
