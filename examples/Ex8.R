# Anderlucci's meat with UFS
# Anderlucci's: https://rdrr.io/cran/RPEClust/f/
# UFS: https://github.com/farhadabedinzadeh/AutoUFSTool

library(R.matlab)
library(diceRplus)
library(RPEClust)
library(future.apply)
library(mclust)
library(clue)

# Original RPGMMClu in Anderlucci's RPEClust
RPGMMClu<-function(x, true.cl=NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, modelNames = NULL, diagonal = FALSE, ensmethod="DWH", seed = 101, verb = FALSE){

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

# Parallel implementation of RPGMMClu
RPGMMClu_parallel2 <- function(x, true.cl=NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, modelNames = NULL, diagonal = FALSE, ensmethod="DWH", seed = 101, verb = FALSE){

  options(future.globals.maxSize = 1 * 1024^3)  # 1 GB

  p <- ncol(x)
  if(is.null(d)) d <- ceiling(c * log(g))
  n <- nrow(x)

  set.seed(seed)
  RPbase <- generateRP(p=ncol(x), d=d, B=B)
  index <- matrix(1:(d*B), d, B)

  plan(cluster, workers = round(availableCores() * 0.8))

  results <- future_lapply(1:B, function(b) {
    A <- as.matrix(RPbase[, index[, b]])
    A.bar <- qr.Q(qr(A), complete=TRUE)[, (d+1):p]
    y <- x %*% A
    y.bar <- x %*% A.bar
    y.star <- cbind(y, y.bar)

    out <- Mclust(y, g, verbose=FALSE, modelNames = modelNames)
    loglik.c <- out$loglik
    cl.m <- out$classification
    bic.c <- out$bic  # BIC for the d "relevant variables", i.e. the d RPs

    # regression model of the (p-d) not relevant variables on the d relevant ones
    X <- cbind(matrix(1, n), y)
    Y <- y.bar
    B <- solve(t(X) %*% X) %*% t(X) %*% Y
    s.ybar.y <- 1 / (n - 1) * t(Y) %*% (diag(n) - X %*% solve(t(X) %*% X, tol=1e-30) %*% t(X)) %*% Y
    if (diagonal) s.ybar.y <- diag(diag(s.ybar.y))

    det.o <- if (diagonal) prod(diag(s.ybar.y)) else det(s.ybar.y)
    det.s.ybar.y <- if ((det.o) < 1e-323) 1e-323 else if (det.o == Inf) 1e+308 else det.o

    loglik.nc <- (-(ncol(y.bar) * n) / 2) * log((2 * pi)) - (n / 2) * log(det.s.ybar.y) - 0.5 * sum(diag(((Y - X %*% B) %*% solve(s.ybar.y, tol=1e-30) %*% t(Y - X %*% B))))

    # BIC for the (p-d) "not relevant variables"
    k <- if (diagonal) ncol(y.bar) * (d + 1) + ncol(y.bar) else ncol(y.bar) * (d + 1) + (ncol(y.bar) * (ncol(y.bar) + 1)) / 2
    bic.nc <- 2 * loglik.nc - k * log(n)

    list(loglik.c=loglik.c, cl.m=cl.m, bic.c=bic.c, bic.nc=bic.nc, Bic=bic.c + bic.nc, Ari=if (!is.null(true.cl)) adjustedRandIndex(out$classification, true.cl) else NA)
  })
  plan(sequential)
  gc()

  # Convert results to vectors
  Bic <- sapply(results, function(x) x$Bic)
  cl.m <- do.call(cbind, lapply(results, function(x) x$cl.m))
  Ari <- sapply(results, function(x) x$Ari)

  # Ensemble
  cl.ens.1 <- data.frame(cl.m[, order(Bic, decreasing=TRUE)[1:B.star]])
  cl.ens.1.2 <- lapply(cl.ens.1, function(x) as.cl_membership(x))
  cl.consensus <- apply(cl_consensus(cl.ens.1.2, method=ensmethod)$.Data, 1, which.max)

  if (!is.null(true.cl)) ari <- adjustedRandIndex(cl.consensus, true.cl)

  names(ari) <- paste0("B.star=", B.star)
  ensemble <- list(ari = ari, label.vec = cl.consensus)
  individual <- list(label.vec = cl.m, ari = Ari, bic = Bic)

  return(list(ensemble = ensemble, individual = individual))
}

Matlab$startServer()
matlab <- Matlab()
print(matlab)

isOpen <- open(matlab)
if (!isOpen)
  print("MATLAB server is not running: waited 30 seconds.")
print(matlab)

setVariable(matlab, X = Meat$x)

evaluate(matlab, "whos")

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

UFS_Methods_Slow <- list(
  "MCFS" = TRUE,
  "LLCFS" = TRUE,
  "FSASL" = TRUE,
  "SOCFS" = TRUE,
  "SOGFS" = TRUE,
  "MCFS" = TRUE
)

UFS_Results <- list()

for (UFS_Method in names(UFS_Methods))  {
  print(UFS_Method)
  send_input <- UFS_Methods[[UFS_Method]]
  cmd <- sprintf('Result = Auto_UFSTool(X, "%s");', UFS_Method)
  evaluate(matlab, "rng(42);")
  # Execute AutoHotkey to send option 1 keystroke (default values) to MATLAB console
  if (send_input) {
    system('cmd /c start /B "C:/Program Files/AutoHotkey/v2/AutoHotkey.exe" "C:/Users/anton/OneDrive/Documentos/send_input.ahk"')
  }
  # Ahora ejecutar la función en MATLAB
  evaluate(matlab, cmd)
  result <- getVariable(matlab, "Result")
  UFS_Results[[UFS_Method]] <- result
}

close(matlab)

str(UFS_Results)



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

print(execution_time)
str(out.clu_p_baseline)


exp_data <- experiment_logger(
  description = "Baseline clustering with RPGMMClu",
  dataset = "Meat",
  ensemble_method = "RPGMMClu",
  ensemble_method_params = list(g = 5, B = 100, B.star = 10),
  UFS_method = "N/A",
  num_features = NA,
  features = list(),
  dim_reduction_method = "N/A",
  dim_reduction_method_params = list(),
  execution_time = as.numeric(execution_time),
  internal_metrics = list(),
  external_metrics = list(ensemble_ari = out.clu_p_baseline$ensemble$ari)
)


save_experiment(exp_data)

loaded_experiments <- load_experiments()
print(loaded_experiments)  # Ver los experimentos almacenados




# Run RPGMMClu baseline in parallel (80% of available cores)
system.time(out.clu_p_baseline_parallel <- RPGMMClu_parallel(Meat$x, Meat$y, g=5, B=B, B.star=B.star, verb=TRUE))


# RPGMMClu with 1000 more important features (UFS=InfFS)

# More relevant features

# Number of features
N <- 2000
# Extract N índices of more important features
top_features <- UFS_Results$InfFS$Result[[1]][1:N]

# Select columns in dataset
Meat$x_filtered <- Meat$x[, top_features]

B <- 1000
B.star=100

# Run RPGMMClu baseline
system.time(out.clu_p_baseline_1000 <- RPGMMClu_parallel(Meat$x, Meat$y, g=5, B=B, B.star=B.star, verb=TRUE))

# Run RPGMMClu on 1000 more important features
system.time(out.clu_p_UFS_1000  <- RPGMMClu_parallel(Meat$x_filtered, Meat$y, g=5, B=B, B.star=B.star, verb=TRUE))


B <- 1500
B.star=150

# Run RPGMMClu baseline
system.time(out.clu_p_baseline_1500 <- RPGMMClu_parallel(Meat$x, Meat$y, g=5, B=B, B.star=B.star, verb=TRUE))

# Run RPGMMClu on 1000 more important features
system.time(out.clu_p_UFS2000_1500  <- RPGMMClu_parallel(Meat$x_filtered, Meat$y, g=5, B=B, B.star=B.star, verb=TRUE))



exp_data <- experiment_logger(
  description = "Clustering evaluation with PCA and InfFS",
  ensemble_method = "Hierarchical Consensus",
  ensemble_method_params = list(threshold = 0.85, iterations = 50),
  UFS_method = "InfFS",
  num_features = 100,
  features = list("feature1", "feature2", "feature3"),
  dim_reduction_method = "PCA",
  dim_reduction_method_params = list(components = 2),
  execution_time = 47.8,  # In seconds
  internal_metrics = list(silhouette = 0.72, dunn = 1.5),
  external_metrics = list(ARI = 0.85, NMI = 0.78)
)

save_experiment(exp_data)  # Store the experiment in an RDS file

loaded_experiments <- load_experiments()  # Load all logged experiments
print(loaded_experiments)  # View stored experiments



