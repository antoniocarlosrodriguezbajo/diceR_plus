library(diceRplus)
data("ALLAML")

UFS_methods <- list(
  "LLCFS"
)

UFS_Results <- runUFS_manual_params(ALLAML$x, UFS_methods)
save(UFS_Results, file = "experiments/UFS_ALLAM_LLCFS.RData")

# load("experiments/UFS_ALLAML.RData")
load("experiments/UFS_ALLAM_LLCFS.RData")

method = "LLCFS"
print(method)
N <- 150
B <- 100
B.star=10
top_features <- UFS_Results$Results[[method]]$Result[[1]][1:N]


execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(ALLAML$x[, top_features],
                                                                    ALLAML$y,
                                                                    g=2,
                                                                    B=B,
                                                                    B.star=B.star,
                                                                    verb=TRUE))["elapsed"]

# Crear la matriz de confusión entre clusters y etiquetas reales
conf_mat <- table(out.clu_baseline$ensemble$label.vec, ALLAML$y)

# Resolver la asignación óptima de etiquetas
optimal_assignment <- solve_LSAP(conf_mat, maximum = TRUE)

# Relabeling: Asignar etiquetas correctas al clustering
pred_labels_reassigned <- optimal_assignment[out.clu_baseline$ensemble$label.vec]
accuracy <- mean(ALLAML$y == pred_labels_reassigned)
accuracy


dice.obj <- dice(
  data = ALLAML$x[, top_features],
  #data = ALLAML$x,
  p.item = 0.85,
  reps = 10,
  nk = 2,
  algorithms = c("nmf"),
  nmf.method = c("brunet"),
  progress = TRUE,
  verbose = FALSE,
  ref.cl = as.integer(as.factor(ALLAML$y))
)

dice.obj$indices$ei

dice.obj$indices$ii

pred_labels <- dice.obj$clusters
conf_stats <- ev_confmat(pred.lab = pred_labels, ref.lab = ALLAML$y)
print(conf_stats)
