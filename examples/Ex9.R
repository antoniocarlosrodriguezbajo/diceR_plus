# Anderlucci's meat with UFS integrated into the package
# Anderlucci's: https://rdrr.io/cran/RPEClust/f/
# UFS: https://github.com/farhadabedinzadeh/AutoUFSTool

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

exp_data <- experiment_logger(
  description = "Baseline clustering with RPGMMClu",
  dataset = "Meat",
  ensemble_method = "RPGMMClu",
  ensemble_method_params = list(g = 5, B = 100, B.star = 10),
  UFS_method = "N/A",
  num_features = NA,
  features = NA,
  dim_reduction_method = "N/A",
  dim_reduction_method_params = NA,
  execution_time = as.numeric(execution_time),
  labels_clustering = out.clu_p_baseline$ensemble$label.vec,
  internal_metrics = NA,
  external_metrics = list(ensemble_ari = out.clu_p_baseline$ensemble$ari)
)

# Ver la estructura del tibble
str(exp_data)



library(tibble)
# Crear un tibble con diferentes tipos de datos
mi_tibble <- tibble(
  ID = 1:5,  # Números enteros
  Nombre = c("Ana", "Luis", "Sofía", "Carlos", "Marta"),  # Caracteres
  Edad = c(28, 32, 25, 40, 22),  # Números enteros
  Activo = c(TRUE, FALSE, TRUE, TRUE, FALSE),  # Valores lógicos
  Puntajes = list(c(80, 90, 85), c(95, 92), c(78, 88, 91, 85), c(85, 89), c(92, 90, 87))  # Lista con vectores de longitud variable
)

# Ver la estructura del tibble
print(mi_tibble)
