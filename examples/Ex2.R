install.packages("apcluster")
install.packages("kernlab")
install.packages("RandPro")
install.packages("poLCA")

install.packages("https://cran.r-project.org/src/contrib/Archive/clusteval/clusteval_0.1.tar.gz", repos=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/RPEClust/RPEClust_0.1.0.tar.gz", repos=NULL, type="source")

library(RPEClust)




data(Meat)

set.seed(seed)
RPbase<-generateRP(p=ncol(Meat$x),d=17,B=1000)


out.clu <- RPGMMClu(Meat$x, Meat$y, g=5, B=100, B.star=10, verb=TRUE)

data <- sim_normal(n = rep(100, 2), p = 100, rho = rep(0.1, 2), delta = 0.5, sigma2 = 1, seed = 106)
out.clu <- RPGMMClu(data$x, data$y, g=2, B=10, B.star=5, verb=TRUE)


library(parallel)
library(foreach)
library(doParallel)
library(mclust)
library(clue) # Para cl_consensus y otras funciones necesarias

library(parallel)
library(foreach)
library(doParallel)
library(mclust) # Para Mclust
library(clue)   # Para cl_consensus y as.cl_membership

library(parallel)
library(foreach)
library(doParallel)
library(mclust) # Para Mclust
library(clue)   # Para cl_consensus y as.cl_membership

RPGMMClu_paral <- function(x, true.cl = NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, modelNames = NULL, diagonal = FALSE, ensmethod = "DWH", seed = 101, verb = FALSE) {
  p <- ncol(x)
  if (is.null(d)) d <- ceiling(c * log(g)) else d <- d
  n <- nrow(x)

  set.seed(seed)
  RPbase <- generateRP(p = ncol(x), d = d, B = B)
  index <- matrix(1:(d * B), d, B)

  # Preparación para la paralelización
  cl <- makeCluster(detectCores() - 1) # Detecta el número de núcleos menos 1
  registerDoParallel(cl)

  # Ejecutar el bucle en paralelo
  results <- foreach(b = 1:B, .combine = rbind, .packages = c("mclust", "clue")) %dopar% {
    # Validar índices
    if (b > ncol(index)) stop("Índice b fuera de los límites")

    # Generar A, y y y.bar
    A <- as.matrix(RPbase[, index[, b]])
    A.bar <- qr.Q(qr(A), complete = TRUE)[, (d + 1):p]
    y <- x %*% A
    y.bar <- x %*% A.bar
    y.star <- cbind(y, y.bar)

    # Modelo GMM
    out <- Mclust(y, g, verbose = FALSE, modelNames = modelNames)
    loglik.c <- out$loglik
    cl.m <- out$classification
    bic.c <- out$bic

    # Modelo de regresión
    X <- cbind(matrix(1, n), y)
    Y <- y.bar
    B <- solve(t(X) %*% X) %*% t(X) %*% Y
    s.ybar.y <- 1 / (n - 1) * t(Y) %*% (diag(n) - X %*% solve(t(X) %*% X, tol = 1e-30) %*% t(X)) %*% Y
    if (diagonal) s.ybar.y <- diag(diag(s.ybar.y))
    m <- ncol(y.bar)
    if (diagonal) det.o <- prod(diag(s.ybar.y)) else det.o <- det(s.ybar.y)
    if ((det.o) < 1e-323) det.s.ybar.y <- 1e-323 else if (det.o == Inf) det.s.ybar.y <- 1e+308 else det.s.ybar.y <- det.o

    loglik.nc <- (-(m * n) / 2) * log((2 * pi)) - (n / 2) * log(det.s.ybar.y) - 0.5 * sum(diag(((Y - X %*% B) %*% solve(s.ybar.y, tol = 1e-30) %*% t(Y - X %*% B))))

    # BIC para las variables no relevantes
    if (diagonal) k <- m * (d + 1) + m else k <- m * (d + 1) + (m * (m + 1)) / 2
    bic.nc <- 2 * loglik.nc - k * log(n)

    # Calcular BIC final
    Bic <- bic.c + bic.nc
    if (!is.null(true.cl)) Ari <- adjustedRandIndex(out$classification, true.cl) else Ari <- NA

    list(loglik.c = loglik.c, cl.m = cl.m, bic.c = bic.c, bic.nc = bic.nc, Bic = Bic, Ari = Ari)
  }

  stopCluster(cl) # Detener el clúster

  # Extraer resultados del foreach
  loglik.c <- sapply(results, `[[`, "loglik.c")
  cl.m <- do.call(cbind, lapply(results, `[[`, "cl.m"))
  bic.c <- sapply(results, `[[`, "bic.c")
  bic.nc <- sapply(results, `[[`, "bic.nc")
  Bic <- sapply(results, `[[`, "Bic")
  Ari <- sapply(results, `[[`, "Ari")

  # Ensemble ####
  cl.ens.1 <- data.frame(cl.m[, order(Bic, decreasing = TRUE)[1:B.star]])
  cl.ens.1.2 <- lapply(cl.ens.1, function(x) as.cl_membership(x))
  cl.consensus <- apply(cl_consensus(cl.ens.1.2, method = ensmethod)$.Data, 1, which.max)
  if (!is.null(true.cl)) ari <- adjustedRandIndex(cl.consensus, true.cl)

  names(ari) <- paste0("B.star=", B.star)
  ensemble <- list(ari = ari, label.vec = cl.consensus)
  individual <- list(label.vec = cl.m, ari = Ari, bic = Bic, bic.GMM = bic.c, bic.reg = bic.nc)
  return(list(ensemble = ensemble, individual = individual))
}

out.clu <- RPGMMClu_paral(Meat$x, Meat$y, g=5, B=100, B.star=10, verb=TRUE)


library(doParallel)
library(foreach)

numCores <- 8
cl <- makeCluster(numCores)
registerDoParallel(cl)
# Ejecutar en paralelo
results <- foreach(i = 1:numCores, .combine = rbind, .packages = "RPEClust") %dopar% {
  out.clu <- RPGMMClu(Meat$x, Meat$y, g=5, B=100, B.star=10, verb=TRUE)
  return(out.clu)
}
stopCluster(cl)


RPGMMClu_orig<-function(x, true.cl=NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, modelNames = NULL, diagonal = FALSE, ensmethod="DWH", seed = 101, verb = FALSE){

  p<-ncol(x)
  if(is.null(d)) d<-ceiling(c*log(g)) else d = d
  n<-nrow(x)

  set.seed(seed)
  RPbase<-generateRP(p=ncol(x),d=d,B=B)
  index=matrix(1:(d*B),d,B)

  Bic<-Ari<-bic.c<-bic.nc<-loglik.c<-det.s.ybar.y<-NULL
  cl.m<-matrix(NA, n, B)

  for (b in 1:B){
    A<-as.matrix(RPbase[,index[,b]])
    A.bar <- qr.Q(qr(A),complete=TRUE)[,(d+1):p]
    y<-x%*%A
    y.bar<-x%*%A.bar
    y.star<-cbind(y,y.bar)
    out<-Mclust(y,g, verbose=FALSE, modelNames = modelNames)
    loglik.c[b]<-out$loglik
    cl.m[,b]<-out$cl
    bic.c[b]<-out$bic # BIC for the d "relevant variables", i.e. the d RPs

    #regression model of the (p-d) not relevant variables on the d relevant ones
    X<-cbind(matrix(1,n),y)
    Y<-y.bar
    B<-solve(t(X)%*%X)%*%t(X)%*%Y
    s.ybar.y<-1/(n-1)*t(Y)%*%(diag(n)-X%*%solve(t(X)%*%X, tol=1e-30)%*%t(X))%*%Y
    if(diagonal) s.ybar.y<-diag(diag(s.ybar.y))
    m<-ncol(y.bar)
    if(diagonal) det.o<-prod(diag(s.ybar.y)) else det.o<-det(s.ybar.y)
    if((det.o)<1e-323) det.s.ybar.y[b]<-1e-323 else if(det.o==Inf) det.s.ybar.y[b]<-1e+308 else det.s.ybar.y[b]<-det.o

    loglik.nc<-(-(m*n)/2)*log((2*pi))-(n/2)*log(det.s.ybar.y[b])-0.5*sum(diag(((Y-X%*%B)%*%solve(s.ybar.y,tol=1e-30)%*%t(Y-X%*%B))))

    # BIC for the (p-d) "not relevant variables"
    if (diagonal) k<-m*(d+1)+m else k<-m*(d+1)+(m*(m+1))/2  # free parameters for the regression model
    bic.nc[b]<-2*loglik.nc-k*log(n)


    Bic[b]<-bic.c[b]+bic.nc[b]
    if(!is.null(true.cl)) Ari[b]=adjustedRandIndex(out$classification,true.cl)

    if(verb) print(b)
  }

  # Ensemble ####

  cl.ens.1<-data.frame(cl.m[,order(Bic,decreasing=TRUE)[1:B.star]])
  cl.ens.1.2<-lapply(cl.ens.1,function(x) as.cl_membership(x))
  cl.consensus<-apply(cl_consensus(cl.ens.1.2,method=ensmethod)$.Data,1,which.max)
  if(!is.null(true.cl)) ari<-adjustedRandIndex(cl.consensus,true.cl)


  names(ari)<-paste0("B.star=",B.star)
  ensemble<-list(ari = ari, label.vec = cl.consensus)
  individual<-list(label.vec = cl.m, ari = Ari, bic = Bic, bic.GMM=bic.c, bic.reg=bic.nc)
  return(list(ensemble = ensemble, individual = individual))
}

out.clu <- RPGMMClu_orig(Meat$x, Meat$y, g=5, B=10, B.star=1, verb=TRUE)

out.clu <- RPGMMClu_paral(Meat$x, Meat$y, g=5, B=10, B.star=1, verb=TRUE)


library(foreach)
library(doParallel)
library(mclust)  # Para Mclust

library(foreach)
library(doParallel)
library(mclust)  # Para Mclust

library(foreach)
library(doParallel)
library(mclust)  # Para Mclust

library(foreach)
library(doParallel)
library(mclust)

library(foreach)
library(doParallel)
library(mclust)

library(foreach)
library(doParallel)
library(mclust)

RPGMMClu_parallel <- function(x, true.cl=NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, modelNames = NULL, diagonal = FALSE, ensmethod="DWH", seed = 101, verb = FALSE) {

  p <- ncol(x)
  if (is.null(d)) d <- ceiling(c * log(g)) else d <- d
  n <- nrow(x)

  set.seed(seed)
  RPbase <- generateRP(p = ncol(x), d = d, B = B)
  index <- matrix(1:(d * B), d, B)

  # Configuración de la paralelización
  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  # Ejecutar el bucle en paralelo
  results <- foreach(b = 1:B, .packages = c("mclust"), .combine = list) %dopar% {

    A <- as.matrix(RPbase[, index[, b]])
    A.bar <- qr.Q(qr(A), complete = TRUE)[, (d + 1):p]
    y <- x %*% A
    y.bar <- x %*% A.bar
    y.star <- cbind(y, y.bar)

    out <- Mclust(y, g, verbose = FALSE, modelNames = modelNames)

    list(loglik.c = out$loglik, bic.c = out$bic, classification = out$classification)
  }

  stopCluster(cl)

  # Convertir resultados en un dataframe
  results_df <- do.call(rbind, results)

  # Reconstrucción segura de `cl.m`
  cl.m <- matrix(unlist(results_df$classification), ncol = B, byrow = TRUE)

  # Seleccionar mejores B.star resultados
  B.star <- min(B.star, ncol(cl.m))
  cl.ens.1 <- data.frame(cl.m[, order(results_df$bic.c, decreasing = TRUE)[1:B.star]])
  cl.ens.1.2 <- lapply(cl.ens.1, function(x) as.cl_membership(x))

  # Evitar errores de dimensión en `cl_consensus`
  cl.consensus <- if (nrow(cl.ens.1) > 0) {
    apply(cl_consensus(cl.ens.1.2, method = ensmethod)$.Data, 1, which.max)
  } else {
    rep(NA, n)
  }

  # Ajustar la comparación con `true.cl`
  ari <- if (!is.null(true.cl) && length(cl.consensus) == length(true.cl)) {
    adjustedRandIndex(cl.consensus, true.cl)
  } else {
    NA
  }

  names(ari) <- paste0("B.star=", B.star)
  ensemble <- list(ari = ari, label.vec = cl.consensus)
  individual <- list(label.vec = cl.m, ari = results_df$classification, bic = results_df$bic.c, bic.GMM = results_df$bic.c, bic.reg = results_df$bic.nc)

  return(list(ensemble = ensemble, individual = individual))
}



out.clu <- RPGMMClu(Meat$x, Meat$y, g=5, B=100, B.star=10, verb=TRUE)


library(mclust)
library(clue)
library(clusteval)

system.time(out.clu <- RPGMMClu(Meat$x, Meat$y, g=5, B=100, B.star=10, verb=TRUE))
system.time(out.clu <- RPGMMClu_parallel(Meat$x, Meat$y, g=5, B=100, B.star=10, verb=TRUE))


library(apcluster)
library(diceR)
library(kernlab)
library(mclust)
library(RandPro)
library(poLCA)

data(meats)

dat <-  meats[, !names(meats) %in% "class"]

ref.cl <- as.integer(meats$class)

meats <- dat

dim(meats)

# Parámetros
n <- nrow(meats)
p <- ncol(meats)
G <- 5
d <- ceiling(10 * log(G) + 1)  # ≈ 17, como en el artículo
d <- 17

# Matriz de proyección aleatoria
set.seed(123)
proj_matrix <- form_matrix(p, d, "gaussian")

# Convierte a matriz numérica
meats_mat <- as.matrix(meats)

# Proyección de los datos
meats_rp <- meats_mat %*% proj_matrix
dim(meats_rp) # 231 x 17

meats_rp_scaled <- prepare_data(meats_rp, scale = TRUE)

result <- dice(meats_rp_scaled, nk = 5, reps = 100,
               cons.funs = "kmodes",
               algorithms = c("hc", "km", "pam", "ap", "sc", "gmm"))

cluster_assignments <- result$clusters

lce_clustering <- cluster_assignments[, "LCE"]
clusters_relabel <- relabel_class(lce_clustering,ref.cl)
adjustedRandIndex(clusters_relabel, ref.cl)

clusters_relabel <- relabel_class(cluster_assignments,ref.cl)
adjustedRandIndex(clusters_relabel, ref.cl)




colMeans(dat)
apply(dat, 2, sd)
boxplot(colMeans(dat), main="Medias por variable")
boxplot(apply(dat, 2, sd), main="SD por variable")



set.seed(1) # Para reproducibilidad

dice.obj <- dice(
  data = dat,
  nk = 5,
  reps = 100, # Puedes aumentar a 1000 si tienes suficiente poder de cómputo
  algorithms = c("hc", "km", "pam", "ap", "sc", "gmm"),
  cons.funs = "kmodes",
  progress = TRUE,
  verbose = TRUE
)

# Para ver los clusters asignados:
clusters <- dice.obj$clusters
table(clusters)

clusters_relabel <- relabel_class(clusters,ref.cl)


adjustedRandIndex(clusters_relabel, ref.cl)


library(parallel)
library(mclust)
# Asegúrate de definir/poner aquí tus funciones auxiliares: generateRP, as.cl_membership, cl_consensus, adjustedRandIndex

RPGMMClu_par <- function(x, true.cl=NULL, g, d = NULL, c = 10, B = 1000, B.star = 100,
                     modelNames = NULL, diagonal = FALSE, ensmethod="DWH", seed = 101, verb = FALSE, n.cores = 2) {

  p <- ncol(x)
  if(is.null(d)) d <- ceiling(c*log(g)) else d = d
  n <- nrow(x)

  set.seed(seed)
  RPbase <- generateRP(p=ncol(x), d=d, B=B)
  index <- matrix(1:(d*B), d, B)

  # Función individual para cada iteración
  process_b <- function(b) {
    A <- as.matrix(RPbase[, index[,b]])
    A.bar <- qr.Q(qr(A), complete=TRUE)[, (d+1):p]
    y <- x %*% A
    y.bar <- x %*% A.bar
    y.star <- cbind(y, y.bar)
    out <- Mclust(y, g, verbose=FALSE, modelNames=modelNames)
    loglik_c <- out$loglik
    cl_b <- out$cl
    bic_c <- out$bic

    X <- cbind(matrix(1, n), y)
    Y <- y.bar
    Bmat <- solve(t(X)%*%X) %*% t(X) %*% Y
    s.ybar.y <- 1/(n-1) * t(Y) %*% (diag(n) - X %*% solve(t(X)%*%X, tol=1e-30) %*% t(X)) %*% Y
    if(diagonal) s.ybar.y <- diag(diag(s.ybar.y))
    m <- ncol(y.bar)
    if(diagonal) det.o <- prod(diag(s.ybar.y)) else det.o <- det(s.ybar.y)
    if((det.o) < 1e-323) det_s_ybar_y <- 1e-323 else if(det.o==Inf) det_s_ybar_y <- 1e+308 else det_s_ybar_y <- det.o
    loglik_nc <- (-(m*n)/2)*log((2*pi)) - (n/2)*log(det_s_ybar_y) - 0.5*sum(diag(((Y-X%*%Bmat)%*%solve(s.ybar.y, tol=1e-30)%*%t(Y-X%*%Bmat))))
    if (diagonal) k <- m*(d+1)+m else k <- m*(d+1)+(m*(m+1))/2
    bic_nc <- 2*loglik_nc - k*log(n)
    Bic <- bic_c + bic_nc
    if(!is.null(true.cl)) Ari <- adjustedRandIndex(out$classification, true.cl) else Ari <- NA

    list(
      loglik_c = loglik_c,
      cl_b = cl_b,
      bic_c = bic_c,
      bic_nc = bic_nc,
      det_s_ybar_y = det_s_ybar_y,
      Bic = Bic,
      Ari = Ari
    )
  }

  # 1. Crea el cluster
  cl <- makeCluster(n.cores)

  # 2. Exporta todos los objetos y funciones necesarias a los workers
  clusterExport(cl, varlist=c("x", "g", "d", "c", "RPbase", "index",
                              "generateRP", "Mclust", "as.cl_membership",
                              "cl_consensus", "adjustedRandIndex",
                              "modelNames", "diagonal", "ensmethod", "true.cl", "seed"),
                envir=environment())

  # 3. Carga los paquetes necesarios en cada worker
  clusterEvalQ(cl, library(mclust))

  # 4. Ejecuta el bucle en paralelo
  results <- parLapply(cl, 1:B, process_b)

  # 5. Detén el cluster
  stopCluster(cl)

  # 6. Recoge resultados
  loglik.c <- sapply(results, `[[`, "loglik_c")
  cl.m <- do.call(cbind, lapply(results, `[[`, "cl_b"))
  bic.c <- sapply(results, `[[`, "bic_c")
  bic.nc <- sapply(results, `[[`, "bic_nc")
  det.s.ybar.y <- sapply(results, `[[`, "det_s_ybar_y")
  Bic <- sapply(results, `[[`, "Bic")
  Ari <- sapply(results, `[[`, "Ari")

  # Ensemble
  cl.ens.1 <- data.frame(cl.m[, order(Bic, decreasing=TRUE)[1:B.star]])
  cl.ens.1.2 <- lapply(cl.ens.1, function(x) as.cl_membership(x))
  cl.consensus <- apply(cl_consensus(cl.ens.1.2, method=ensmethod)$.Data, 1, which.max)
  if(!is.null(true.cl)) ari <- adjustedRandIndex(cl.consensus, true.cl)
  names(ari) <- paste0("B.star=",B.star)
  ensemble <- list(ari = ari, label.vec = cl.consensus)
  individual <- list(label.vec = cl.m, ari = Ari, bic = Bic, bic.GMM = bic.c, bic.reg = bic.nc)
  return(list(ensemble = ensemble, individual = individual))
}

system.time(out.clu <- RPGMMClu_parallel(Meat$x, Meat$y, g=5, B=100, B.star=10, verb=TRUE))
