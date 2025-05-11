# Monti's leukemia

install.packages("kohonen")


library(kohonen)
library(diceRplus)
library(mclust)


normalize_rows_and_cols <- function(M) {
  M1 <- t(scale(t(M), center = TRUE, scale = TRUE)) # filas
  M2 <- scale(M1, center = TRUE, scale = TRUE)      # columnas
  return(M2)
}


ks <- c(3)

data(leukemia)

dat <-  leukemia[, !names(leukemia) %in% "class"]


ref.cl <- as.integer(leukemia$class)

top_n_genes <- 999
vars <- apply(dat, 2, var)
selected_genes <- order(vars, decreasing = TRUE)[1:top_n_genes]
dat_selected <- dat[, selected_genes]

dat_selected <- normalize_rows_and_cols(dat_selected)

set.seed(123)
resultado <- dice(
  data = dat_selected,
  nk = 3,                   # probar varios valores de K
  reps = 500,                 # número de resamples recomendado
  algorithms = "hc",          # solo hierarchical clustering
  hc.method = "average",      # linkage promedio, como en el paper
  cons.funs = c("majority"),    # método de combinación (kmodes ≈ majority voting)
  evaluate = TRUE,            # calcula índices de validez (incluye PAC, Rand, etc.)
  ref.cl = ref.cl,
  # plot = TRUE,                # genera gráficas de CDF, PAC, etc.
  progress = TRUE
)

resultado$indices$ei

# Supón que leukemia es un data.frame o matriz
apply(dat_norm, 1, mean)    # medias de filas
apply(dat_norm, 2, mean)    # medias de columnas
apply(dat_norm, 1, sd)      # sd de filas
apply(dat_norm, 2, sd)      # sd de columnas
