install.packages("kohonen")


library(kohonen)
library(diceRplus)
library(mclust)
ks <- c(3)

data(leukemia)

dat <-  leukemia[, !names(leukemia) %in% "class"]
ref.cl <- as.integer(leukemia$class)

top_n_genes <- 999
vars <- apply(dat, 2, var)
selected_genes <- order(vars, decreasing = TRUE)[1:top_n_genes]
dat_selected <- dat[, selected_genes]

set.seed(123)
resultado <- dice(
  data = dat_selected,
  nk = 3,                   # probar varios valores de K
  reps = 500,                 # número de resamples recomendado
  algorithms = "hc",          # solo hierarchical clustering
  hc.method = "average",      # linkage promedio, como en el paper
  cons.funs = c("kmodes"),    # método de combinación (kmodes ≈ majority voting)
  evaluate = TRUE,            # calcula índices de validez (incluye PAC, Rand, etc.)
  ref.cl = ref.cl,
  # plot = TRUE,                # genera gráficas de CDF, PAC, etc.
  progress = TRUE
)

resultado$indices$ei

