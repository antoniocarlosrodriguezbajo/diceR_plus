data(meats)

str(meats$class)

dat <-  meats[, !names(meats) %in% "class"]

ref.cl <- as.integer(meats$class)

set.seed(123)
dice.obj <- dice(
  dat,
  nk = 5,
  reps = 15,
  algorithms = c("hc","km"),
  cons.funs = c("majority", "CSPA"),
  progress = TRUE,
  verbose = TRUE
)
k_optimo <- dice.obj$indices$k


# Supongamos que quieres relabelar los clusters obtenidos por majority voting:
clusters <- dice.obj$clusters # Esto es un vector con los labels de cluster

# Relabel para que los labels correspondan lo más posible a los de ref.cl:
clusters_relabel <- relabel_class(clusters[, 1], ref.cl)

# Ahora clusters_relabel tiene los mismos valores de etiqueta que ref.cl (en la medida de lo posible)
table(ref.cl, clusters_relabel)

ev_nmi(clusters_relabel, ref.cl)


library(mclust)
adjustedRandIndex(clusters_relabel, ref.cl)

adjustedRandIndex(clusters[, 2], ref.cl)


# Acceder a la lista de índices internos:
indices <- dice.obj$indices

# Ver todas las métricas internas para el k óptimo
indices$ii[[as.character(k_optimo)]]

# Ver todas las métricas externas para el k óptimo
indices$ei[[as.character(k_optimo)]]


names(dice.obj$indices)

# PAC del ensemble final:
dice.obj$indices$pac

# Asignación final de clusters según el ensemble:
dice.obj$clusters

dice.obj$clusters
table(dice.obj$clusters)
str(dice.obj, max.level = 2)

# Acceder a la lista de métricas internas para k=5 (en tu ejemplo)
dice.obj$indices$ii$`5`

dice.obj$indices$ii


dice.obj <- dice(
  dat,
  nk = 5,
  reps = 15,
  algorithms = c("hc","km"),
  cons.funs = c("majority", "CSPA"),
  progress = TRUE,
  verbose = TRUE
)
dice.obj$indices$ii$`5`




















data(meats)

dat <-  meats[, !names(meats) %in% "class"]

set.seed(123)
dice.obj <- dice(
  dat,
  nk = 4:6,
  reps = 15,
  algorithms = c("hc","km"),
  cons.funs = c("majority", "CSPA"),
  progress = TRUE,
  verbose = TRUE
)

# Valor óptimo de k elegido por el ensemble:
dice.obj$k

# PAC del ensemble final:
dice.obj$indices$pac

# Saber cuál es el k óptimo elegido por el ensemble:
k_optimo <- dice.obj$indices$k

# Acceder a la lista de índices internos:
indices <- dice.obj$indices

# Ver todas las métricas internas para el k óptimo
indices$ii[[as.character(k_optimo)]]


set.seed(123)
cc <- consensus_cluster(
  data = meats_features,
  nk = 5,                  # Número de clusters deseado
  reps = 10,               # Número de repeticiones/submuestras
  algorithms = c("km", "pam"),  # K-means y Partition Around Medoids
  progress = TRUE
)
# Obtener las clases de consenso usando la función k-modes
clases_ensamble_kmodes <- k_modes(cc[, , 1, 1, drop = FALSE])
table(clases_ensamble_kmodes)

# Obtener las clases de consenso usando la función majority_voting
clases_ensamble_majority <- majority_voting(cc[, , 1, 1, drop = FALSE])
table(clases_ensamble_kmodes)

score_kmodes <- compactness(meats_features, clases_ensamble_kmodes)
score_majority <- compactness(meats_features, clases_ensamble_majority)
print(score_kmodes)
print(score_majority)

# Evaluar ambos ensambles con consensus_evaluate
eval_kmodes <- consensus_evaluate(
  meats_features,
  cons.cl = matrix(clases_ensamble_kmodes, ncol = 1)
)
eval_majority <- consensus_evaluate(
  meats_features,
  cons.cl = matrix(clases_ensamble_majority, ncol = 1)
)

# Ver resultados (índices internos)
str(eval_kmodes\$ii, max.level = 2)
str(eval_majority\$ii, max.level = 2)



cc <- consensus_cluster(meats_features, nk = 5, reps = 10, algorithms =c("km", "pam"), progress = FALSE)
clases_ensamble_kmodes <- k_modes(cc[, , 1, 1, drop = FALSE])
eval_kmodes <- consensus_evaluate(
  meats_features,
  cc, # <--- ESTE ES EL FALTANTE
  cons.cl = matrix(clases_ensamble_kmodes, ncol = 1)
)
eval_kmodes <- consensus_evaluate(
  meats_features,
  cc,  # objeto de consensus_cluster
  cons.cl = matrix(clases_ensamble_kmodes, ncol = 1)
)

# Si solo hay un k:
eval_kmodes$ii[[1]]


cc <- consensus_cluster(meats_features, nk = 5, reps = 10, algorithms =c("km", "pam"), progress = FALSE)
clases_ensamble_kmodes <- k_modes(cc[, , 1, 1, drop = FALSE])
eval_kmodes <- consensus_evaluate(
  meats_features,
  cc, # <--- ESTE ES EL FALTANTE
  cons.cl = matrix(clases_ensamble_kmodes, ncol = 1)
)
eval_kmodes <- consensus_evaluate(
  meats_features,
  cc,  # objeto de consensus_cluster
  cons.cl = matrix(clases_ensamble_kmodes, ncol = 1)
)

# Si solo hay un k:
eval_kmodes$ii[[1]]


cc <- consensus_cluster(meats_features, nk = 5, reps = 10, algorithms =c("km", "pam"), progress = FALSE)
clases_ensamble_kmodes <- k_modes(cc[, , 1, 1, drop = FALSE])
eval_kmodes <- consensus_evaluate(
  meats_features,
  cc, # <--- ESTE ES EL FALTANTE
  cons.cl = matrix(clases_ensamble_kmodes, ncol = 1)
)
eval_kmodes <- consensus_evaluate(
  meats_features,
  cc,  # objeto de consensus_cluster
  cons.cl = matrix(clases_ensamble_kmodes, ncol = 1)
)

cc <- consensus_cluster(meats_features, nk = 5, reps = 10, algorithms =c("km", "pam"), progress = FALSE)
clases_ensamble_kmodes <- k_modes(cc[, , 1, 1, drop = FALSE])
eval_kmodes <- consensus_evaluate(
  meats_features,
  cc, # <--- ESTE ES EL FALTANTE
  cons.cl = matrix(clases_ensamble_kmodes, ncol = 1)
)
eval_kmodes <- consensus_evaluate(
  meats_features,
  cc,  # objeto de consensus_cluster
  cons.cl = matrix(clases_ensamble_kmodes, ncol = 1)
)

dat <- meats_features  # Asumiendo que ya está en tu entorno
# Vector de algoritmos de clustering a usar (puedes ajustar)
algorithms <- c("pam", "hc", "km")  # PAM, Hierarchical, K-means

# Lista de métodos ensemble a comparar (deben ser los nombres de las funciones)
ensemble_methods <- list(
  kmodes   = function(E, k) k_modes(E),
  majority = function(E, k) majority_voting(E),
  cspa     = function(E, k) CSPA(E, k = k),
  lce      = function(E, k) LCE(E, k = k, sim.mat = "cts"),
  lca      = function(E, k) LCA(E)
)

set.seed(123)
cc <- consensus_cluster(dat, nk = 5, reps = 6, algorithms = algorithms, progress = TRUE)

# Extraemos la matriz E: todas las muestras x repeticiones x todos los algoritmos x k=4
E <- cc[, , , 1, drop = FALSE]  # El último índice es para k=4

# Inicializamos listas para almacenar resultados
clustering_results <- list()
evaluation_results <- list()

# Extraemos la matriz E: todas las muestras x repeticiones x todos los algoritmos x k=4
E <- cc[, , , 1, drop = FALSE]  # El último índice es para k=4

# Inicializamos listas para almacenar resultados
clustering_results <- list()
evaluation_results <- list()

for (method in names(ensemble_methods)) {
  # Calculamos el ensemble
  if (method %in% c("cspa", "lce")) {
    cl <- ensemble_methods[[method]](E, k = 5)
  } else {
    cl <- ensemble_methods[[method]](E, k = 5)
  }
  # Aseguramos tipo correcto
  cl <- as.integer(cl)
  clustering_results[[method]] <- cl
  evaluation_results[[method]] <- consensus_evaluate(dat, cc, cons.cl = matrix(cl, ncol = 1))
}



ref.cl <- strsplit(rownames(dat), "_") |>
  purrr::map_chr(2) |>
  factor() |>
  as.integer()
dice.obj <- dice(dat, nk = 4, reps = 5, algorithms = "hc", cons.funs =
                   "kmodes", ref.cl = ref.cl, progress = FALSE, verbose = FALSE)
str(dice.obj, max.level = 2)


data(meats)

dat <-  meats[, !names(meats) %in% "class"]
dice.obj <- dice(dat, nk = 5, reps = 5, algorithms = "hc",
                 cons.funs = "kmodes", progress = TRUE, verbose = TRUE)

dice.obj$clusters
table(dice.obj$clusters)
str(dice.obj, max.level = 2)

# Acceder a la lista de métricas internas para k=5 (en tu ejemplo)
dice.obj$indices$ii$`5`

dice.obj$indices$ii


dice.obj <- dice(
  dat,
  nk = 5,
  reps = 15,
  algorithms = c("hc","km"),
  cons.funs = c("majority", "CSPA"),
  progress = TRUE,
  verbose = TRUE
)
dice.obj$indices$ii$`5`
