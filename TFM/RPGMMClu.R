#' @title Random Projection Ensemble Clustering Algorithm
#'
#' @description This function allows to run the RPEClu algorithm.
#'
#' @param x A numeric high-dimensional matrix where rows correspond to observations and columns correspond to variables.
#' @param true.cl A vector containing the true cluster membership. If supplied, the Adjusted Rand Index (ARI) of the predicted clustering is also returned. By default is set to NULL.
#' @param g The number of clusters.
#' @param d The dimension of the projected space. If is \code{NULL} (default option), then \eqn{d = \lceil c* log(g) \rceil.}{d = ceil(c* log(g)).}
#' @param c The constant which governs the dimension of the projected space if \eqn{d} is not provided. The default is set to 10.
#' @param B The number of generated projections; the default is 1000.
#' @param B.star The number of \emph{base} models to retain in the final ensemble; the default is 100.
#' @param modelNames A vector of character strings indicating the models to be fitted in the EM phase of the GMM. The help file for \link[mclust]{Mclust} describes the available options.
#' @param seed A single value indicating the initializing seed for random generations.
#' @param diagonal A logical value indicating whether the conditional covariate matrix has a restricted form, i.e. it is diagonal. The default is FALSE.
#' @param ensmethod A character string specifying the method for computing the clustering consensus. See the \link[clue]{cl_consensus} help file for available options. The default is \code{DWH}.
#' @param verb A logical controlling if the progress number of projections is displayed during the fitting procedure. By default is \code{FALSE}.
#' @return The output components are as follows:
#' \item{\code{ensemble}}{A list including:
#' \describe{
#'  \item{\code{label.vec}}{The vector of labels predicted by the ensemble of size \code{B.star}.}
#'  \item{\code{ari}}{The corresponding ARI (if \code{true.cl} is not \code{NULL}).}}}
#' \item{individual}{A list including:
#'  \describe{\item{\code{label.vec}}{The vectors of labels predicted by each \emph{base} model.}
#'   \item{\code{ari}}{The corresponding ARI (if \code{true.cl} is not \code{NULL}).}
#'   \item{\code{bic}}{The BIC associated to each \emph{base} model.}
#'   \item{\code{bic.GMM}}{The BIC associated to the Gaussian mixture fitted on each projected data.}
#'   \item{\code{bic.reg}}{The BIC for the linear regression of the \eqn{(p-d)} last columns of \eqn{Y*} on the first \eqn{d} ones.}}}
#' @references Anderlucci, Fortunato, Montanari (2019) <arXiv:1909.10832>
#' @examples
#' \donttest{
#' data(Meat)
#' out.clu <- RPGMMClu(Meat$x, Meat$y, g=5, B=1000, B.star=100, verb=TRUE)}
#'
#' data <- sim_normal(n = rep(100, 2), p = 100, rho = rep(0.1, 2), delta = 0.5, sigma2 = 1, seed = 106)
#' out.clu <- RPGMMClu(data$x, data$y, g=2, B=10, B.star=5, verb=TRUE)
#'
#' @export

# Anderlucci's: https://rdrr.io/cran/RPEClust/f/
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
    # y.star<-cbind(y,y.bar)
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

  ari <- NULL
  if(!is.null(true.cl)) {
    ari<-adjustedRandIndex(cl.consensus,true.cl)
    names(ari)<-paste0("B.star=",B.star)
  }
  ensemble<-list(ari = ari, label.vec = cl.consensus)
  individual<-list(label.vec = cl.m, ari = Ari, bic = Bic, bic.GMM=bic.c, bic.reg=bic.nc)
  return(list(ensemble = ensemble, individual = individual))
}

#' @title Random Projection Ensemble Clustering Algorithm (Windows parallel implementation)
#'
#' @description This function allows to run the RPEClu algorithm in parallel.
#'
#' @param x A numeric high-dimensional matrix where rows correspond to observations and columns correspond to variables.
#' @param true.cl A vector containing the true cluster membership. If supplied, the Adjusted Rand Index (ARI) of the predicted clustering is also returned. By default is set to NULL.
#' @param g The number of clusters.
#' @param d The dimension of the projected space. If is \code{NULL} (default option), then \eqn{d = \lceil c* log(g) \rceil.}{d = ceil(c* log(g)).}
#' @param c The constant which governs the dimension of the projected space if \eqn{d} is not provided. The default is set to 10.
#' @param B The number of generated projections; the default is 1000.
#' @param B.star The number of \emph{base} models to retain in the final ensemble; the default is 100.
#' @param modelNames A vector of character strings indicating the models to be fitted in the EM phase of the GMM. The help file for \link[mclust]{Mclust} describes the available options.
#' @param seed A single value indicating the initializing seed for random generations.
#' @param diagonal A logical value indicating whether the conditional covariate matrix has a restricted form, i.e. it is diagonal. The default is FALSE.
#' @param ensmethod A character string specifying the method for computing the clustering consensus. See the \link[clue]{cl_consensus} help file for available options. The default is \code{DWH}.
#' @param verb A logical controlling if the progress number of projections is displayed during the fitting procedure. By default is \code{FALSE}.
#' @return The output components are as follows:
#' \item{\code{ensemble}}{A list including:
#' \describe{
#'  \item{\code{label.vec}}{The vector of labels predicted by the ensemble of size \code{B.star}.}
#'  \item{\code{ari}}{The corresponding ARI (if \code{true.cl} is not \code{NULL}).}}}
#' \item{individual}{A list including:
#'  \describe{\item{\code{label.vec}}{The vectors of labels predicted by each \emph{base} model.}
#'   \item{\code{ari}}{The corresponding ARI (if \code{true.cl} is not \code{NULL}).}
#'   \item{\code{bic}}{The BIC associated to each \emph{base} model.}
#'   \item{\code{bic.GMM}}{The BIC associated to the Gaussian mixture fitted on each projected data.}
#'   \item{\code{bic.reg}}{The BIC for the linear regression of the \eqn{(p-d)} last columns of \eqn{Y*} on the first \eqn{d} ones.}}}
#' @references Anderlucci, Fortunato, Montanari (2019) <arXiv:1909.10832>
#' @examples
#' \donttest{
#' data(Meat)
#' out.clu <- RPGMMClu_parallel(Meat$x, Meat$y, g=5, B=1000, B.star=100, verb=TRUE)}
#'
#' data <- sim_normal(n = rep(100, 2), p = 100, rho = rep(0.1, 2), delta = 0.5, sigma2 = 1, seed = 106)
#' out.clu <- RPGMMClu(data$x, data$y, g=2, B=10, B.star=5, verb=TRUE)
#'
#' @export

# Parallel implementation of RPGMMClu
RPGMMClu_parallel <- function(x, true.cl=NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, modelNames = NULL, diagonal = FALSE, ensmethod="DWH", seed = 101, verb = FALSE){

  options(future.globals.maxSize = 1 * 16 * 1024^3)  # 1 GB

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
    # y.star <- cbind(y, y.bar)

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

  ari <- NULL
  if (!is.null(true.cl)) {
    ari <- adjustedRandIndex(cl.consensus, true.cl)
    names(ari) <- paste0("B.star=", B.star)
  }
  ensemble <- list(ari = ari, label.vec = cl.consensus)
  individual <- list(label.vec = cl.m, ari = Ari, bic = Bic)

  return(list(ensemble = ensemble, individual = individual))
}

#' @title Random Projection Ensemble Clustering Algorithm (Windows parallel implementation)
#'
#' @description This function allows to run the RPEClu algorithm in parallel.
#'
#' @param x A numeric high-dimensional matrix where rows correspond to observations and columns correspond to variables.
#' @param true.cl A vector containing the true cluster membership. If supplied, the Adjusted Rand Index (ARI) of the predicted clustering is also returned. By default is set to NULL.
#' @param g The number of clusters.
#' @param d The dimension of the projected space. If is \code{NULL} (default option), then \eqn{d = \lceil c* log(g) \rceil.}{d = ceil(c* log(g)).}
#' @param c The constant which governs the dimension of the projected space if \eqn{d} is not provided. The default is set to 10.
#' @param B The number of generated projections; the default is 1000.
#' @param B.star The number of \emph{base} models to retain in the final ensemble; the default is 100.
#' @param modelNames A vector of character strings indicating the models to be fitted in the EM phase of the GMM. The help file for \link[mclust]{Mclust} describes the available options.
#' @param seed A single value indicating the initializing seed for random generations.
#' @param diagonal A logical value indicating whether the conditional covariate matrix has a restricted form, i.e. it is diagonal. The default is FALSE.
#' @param ensmethod A character string specifying the method for computing the clustering consensus. See the \link[clue]{cl_consensus} help file for available options. The default is \code{DWH}.
#' @param verb A logical controlling if the progress number of projections is displayed during the fitting procedure. By default is \code{FALSE}.
#' @return The output components are as follows:
#' \item{\code{ensemble}}{A list including:
#' \describe{
#'  \item{\code{label.vec}}{The vector of labels predicted by the ensemble of size \code{B.star}.}
#'  \item{\code{ari}}{The corresponding ARI (if \code{true.cl} is not \code{NULL}).}}}
#' \item{individual}{A list including:
#'  \describe{\item{\code{label.vec}}{The vectors of labels predicted by each \emph{base} model.}
#'   \item{\code{ari}}{The corresponding ARI (if \code{true.cl} is not \code{NULL}).}
#'   \item{\code{bic}}{The BIC associated to each \emph{base} model.}
#'   \item{\code{bic.GMM}}{The BIC associated to the Gaussian mixture fitted on each projected data.}
#'   \item{\code{bic.reg}}{The BIC for the linear regression of the \eqn{(p-d)} last columns of \eqn{Y*} on the first \eqn{d} ones.}}}
#' @references Anderlucci, Fortunato, Montanari (2019) <arXiv:1909.10832>
#' @examples
#' \donttest{
#' data(Meat)
#' out.clu <- RPGMMClu_parallel(Meat$x, Meat$y, g=5, B=1000, B.star=100, verb=TRUE)}
#'
#' data <- sim_normal(n = rep(100, 2), p = 100, rho = rep(0.1, 2), delta = 0.5, sigma2 = 1, seed = 106)
#' out.clu <- RPGMMClu(data$x, data$y, g=2, B=10, B.star=5, verb=TRUE)
#'
#' @export

# Parallel implementation of RPGMMClu
RPGMMClu_parallel_GPU <- function(x, true.cl=NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, modelNames = NULL, diagonal = FALSE, ensmethod="DWH", seed = 101, verb = FALSE){

  options(future.globals.maxSize = 1 * 16 * 1024^3)  # 1 GB

  p <- ncol(x)
  if(is.null(d)) d <- ceiling(c * log(g))
  n <- nrow(x)

  set.seed(seed)
  RPbase <- generateRP(p=ncol(x), d=d, B=B)
  index <- matrix(1:(d*B), d, B)

  plan(cluster, workers = round(availableCores() * 0.8))
  plan(cluster, workers = 1)

  results <- future_lapply(1:B, function(b, future.seed=TRUE) {

    x_torch <- torch_tensor(x, device = "cuda", dtype = torch_float32())
    A <- as.matrix(RPbase[, index[, b]])
    A.bar <- qr.Q(qr(A), complete=TRUE)[, (d+1):p]
    #y <- x %*% A
    #y.bar <- x %*% A.bar
    # y.star <- cbind(y, y.bar)

    # --- Inicio de la parte GPU con torch ---

    # Convertir las matrices de R a tensores de torch
    # Especifica device = "cuda" y dtype = torch_float() para que se creen en la GPU con float32
    A_torch <- torch_tensor(A, device = "cuda", dtype = torch_float32())
    A_bar_torch <- torch_tensor(A.bar, device = "cuda", dtype = torch_float32())

    # Realizar las multiplicaciones matriciales en la GPU
    y_torch <- torch_mm(x_torch, A_torch)
    y.bar_torch <- torch_mm(x_torch, A_bar_torch)

    # --- Fin de la parte GPU con torch ---

    # Convertir los tensores de torch de vuelta a matrices de R (en la CPU)
    # Esto es necesario para usar Mclust y las otras funciones de R
    y_torch <- y_torch$cpu()
    y <- as.matrix(y_torch)

    y.bar_torch <- y.bar_torch$cpu()
    y.bar <- as.matrix(y.bar_torch)

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

  ari <- NULL
  if (!is.null(true.cl)) {
    ari <- adjustedRandIndex(cl.consensus, true.cl)
    names(ari) <- paste0("B.star=", B.star)
  }
  ensemble <- list(ari = ari, label.vec = cl.consensus)
  individual <- list(label.vec = cl.m, ari = Ari, bic = Bic)

  return(list(ensemble = ensemble, individual = individual))
}

#' @export
#'
# Parallel implementation of RPGMMClu (Block Version)
RPGMMClu_parallel_block <- function(x, true.cl=NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, modelNames = NULL, diagonal = FALSE, ensmethod="DWH", seed = 101, verb = FALSE, blockSize = 10){

  options(future.globals.maxSize = 1 * 16 * 1024^3)  # 1 GB

  p <- ncol(x)
  if(is.null(d)) d <- ceiling(c * log(g))
  n <- nrow(x)

  set.seed(seed)
  RPbase <- generateRP(p=ncol(x), d=d, B=B)

  # Divide the B iterations into blocks
  numBlocks <- ceiling(B / blockSize)
  blockIndices <- split(1:B, rep(1:numBlocks, each = blockSize, length.out = B))

  plan(cluster, workers = round(availableCores() * 0.8))

  results_list <- future_lapply(blockIndices, function(indices_in_block) {
    block_results <- list()
    # Extract and combine projections for the entire block
    A_block <- do.call(cbind, lapply(indices_in_block, function(b) as.matrix(RPbase[, ((b-1)*d + 1):(b*d)])))

    # Perform matrix multiplication for the entire block
    y_block <- x %*% A_block

    # Now, iterate through the block to process each projection set
    for (i in 1:length(indices_in_block)) {
      b <- indices_in_block[i]
      # Extract the projected data for the current b
      y <- y_block[, ((i-1)*d + 1):(i*d)]

      # Need to generate A.bar for each b individually or find a way to batch this as well
      A <- as.matrix(RPbase[, ((b-1)*d + 1):(b*d)]) # Re-extract A to get A.bar
      A.bar <- qr.Q(qr(A), complete=TRUE)[, (d+1):p]
      y.bar <- x %*% A.bar # This is still done individually per b

      out <- Mclust(y, g, verbose=FALSE, modelNames = modelNames)
      loglik.c <- out$loglik
      cl.m <- out$classification
      bic.c <- out$bic

      X <- cbind(matrix(1, n), y)
      Y <- y.bar
      B_reg <- solve(t(X) %*% X) %*% t(X) %*% Y
      s.ybar.y <- 1 / (n - 1) * t(Y) %*% (diag(n) - X %*% solve(t(X) %*% X, tol=1e-30) %*% t(X)) %*% Y
      if (diagonal) s.ybar.y <- diag(diag(s.ybar.y))

      det.o <- if (diagonal) prod(diag(s.ybar.y)) else det(s.ybar.y)
      det.s.ybar.y <- if ((det.o) < 1e-323) 1e-323 else if (det.o == Inf) 1e+308 else det.o

      loglik.nc <- (-(ncol(y.bar) * n) / 2) * log((2 * pi)) - (n / 2) * log(det.s.ybar.y) - 0.5 * sum(diag(((Y - X %*% B_reg) %*% solve(s.ybar.y, tol=1e-30) %*% t(Y - X %*% B_reg))))

      k <- if (diagonal) ncol(y.bar) * (d + 1) + ncol(y.bar) else ncol(y.bar) * (d + 1) + (ncol(y.bar) * (ncol(y.bar) + 1)) / 2
      bic.nc <- 2 * loglik.nc - k * log(n)

      block_results[[as.character(b)]] <- list(loglik.c=loglik.c, cl.m=cl.m, bic.c=bic.c, bic.nc=bic.nc, Bic=bic.c + bic.nc, Ari=if (!is.null(true.cl)) adjustedRandIndex(out$classification, true.cl) else NA)
    }
    block_results
  })

  plan(sequential)
  gc()

  # Combine results from blocks
  all_results <- unlist(results_list, recursive = FALSE)

  Bic <- sapply(all_results, function(x) x$Bic)
  cl.m <- do.call(cbind, lapply(all_results, function(x) x$cl.m))
  Ari <- sapply(all_results, function(x) x$Ari)

  # Ensemble
  cl.ens.1 <- data.frame(cl.m[, order(Bic, decreasing=TRUE)[1:B.star]])
  cl.ens.1.2 <- lapply(cl.ens.1, function(x) as.cl_membership(x))
  cl.consensus <- apply(cl_consensus(cl.ens.1.2, method=ensmethod)$.Data, 1, which.max)

  ari <- NULL
  if (!is.null(true.cl)) {
    ari <- adjustedRandIndex(cl.consensus, true.cl)
    names(ari) <- paste0("B.star=", B.star)
  }
  ensemble <- list(ari = ari, label.vec = cl.consensus)
  individual <- list(label.vec = cl.m, ari = Ari, bic = Bic)

  return(list(ensemble = ensemble, individual = individual))
}

#' @export

RPGMMClu_torch <- function(x, true.cl = NULL, g, d = NULL, c = 10, B = 1000, B.star = 100,
                           modelNames = NULL, diagonal = FALSE, ensmethod = "DWH",
                           seed = 101, verb = FALSE, device = "cuda") {
  p <- ncol(x)
  if (is.null(d)) d <- ceiling(c * log(g))
  n <- nrow(x)

  set.seed(seed)
  RPbase <- generateRP(p = p, d = d, B = B)
  index <- matrix(1:(d * B), d, B)

  Bic <- Ari <- bic.c <- bic.nc <- loglik.c <- det.s.ybar.y <- NULL
  cl.m <- matrix(NA, n, B)

  # Convierte de entrada a tensor en GPU
  x_torch <- torch_tensor(as.matrix(x), device = device)

  for (b in 1:B) {
    # 1. Proyección aleatoria en torch
    A <- torch_tensor(as.matrix(RPbase[, index[, b]]), device = device) # p x d

    # 2. QR (torch): Q es p x d, pero qr completo da p x p
    qr_res <- linalg_qr(A)
    Q <- qr_res[[1]] # p x d

    # --- Chequeo de dimensiones ---
    if ((d + 1) > p) {
      warning(sprintf("Iteración %d: d+1 > p, no hay columnas irrelevantes para proyectar. Se omite.", b))
      next
    }

    # 3. Extrae las columnas irrelevantes (A_bar)
    A_bar <- Q[, (d + 1):p]

    # 4. Proyección relevante e irrelevante
    y_torch <- torch_matmul(x_torch, A)
    y_bar_torch <- torch_matmul(x_torch, A_bar)

    # Chequeo y_bar vacío
    if (y_bar_torch$size(2) == 0) {
      warning(sprintf("Iteración %d: y_bar tiene 0 columnas. Se omite.", b))
      next
    }

    # --- Clustering, pasando a R para Mclust ---
    y <- as.matrix(as_array(y_torch))
    out <- Mclust(y, g, verbose = FALSE, modelNames = modelNames)
    loglik.c[b] <- out$loglik
    cl.m[, b] <- out$cl
    bic.c[b] <- out$bic

    # --- Regresión torch ---
    X_torch <- torch_cat(list(torch_ones(n, 1, device = device), y_torch), dim = 2)
    Y_torch <- y_bar_torch

    XtX <- torch_matmul(X_torch$t(), X_torch)
    # Inversa robusta
    inv_XtX <- tryCatch({
      linalg_inv(XtX)
    }, error = function(e) {
      warning(sprintf("Iteración %d: XtX singular. Se omite.", b))
      next
    })

    XtY <- torch_matmul(X_torch$t(), Y_torch)
    B_torch <- torch_matmul(inv_XtX, XtY)

    diag_n <- torch_eye(n, device = device)
    X_inv_XtX_Xt <- torch_matmul(torch_matmul(X_torch, inv_XtX), X_torch$t())
    s.ybar.y <- (1 / (n - 1)) * torch_matmul(
      Y_torch$t(), torch_matmul(diag_n - X_inv_XtX_Xt, Y_torch)
    )

    if (diagonal) s.ybar.y <- torch_diag(torch_diag(s.ybar.y))

    # Chequeo de matriz de covarianza
    s.ybar.y_r <- as_array(s.ybar.y)
    if (any(dim(s.ybar.y_r) == 0)) {
      warning(sprintf("Iteración %d: s.ybar.y_r es vacía. Se omite.", b))
      next
    }
    if (nrow(s.ybar.y_r) != ncol(s.ybar.y_r)) {
      warning(sprintf("Iteración %d: s.ybar.y_r no es cuadrada. Se omite.", b))
      next
    }

    # Inversa (torch)
    inv_s.ybar.y <- tryCatch({
      linalg_inv(s.ybar.y)
    }, error = function(e) {
      warning(sprintf("Iteración %d: s.ybar.y singular. Se omite.", b))
      next
    })

    # Determinante (torch)
    m <- y_bar_torch$size(2)
    # Si diagonal, producto de diagonal, si no, determinante
    if (diagonal) {
      det.o <- prod(as_array(torch_diag(s.ybar.y)))
    } else {
      det.o <- as_array(linalg_det(s.ybar.y))
    }
    if (det.o < 1e-323) det.s.ybar.y[b] <- 1e-323
    else if (det.o == Inf) det.s.ybar.y[b] <- 1e+308
    else det.s.ybar.y[b] <- det.o

    # --- Logverosimilitud nocluster ---
    Y <- as_array(Y_torch)
    X <- as_array(X_torch)
    B <- as_array(B_torch)
    s.ybar.y_r <- as_array(s.ybar.y)

    loglik.nc <- (-(m * n) / 2) * log((2 * pi)) -
      (n / 2) * log(det.s.ybar.y[b]) -
      0.5 * sum(diag(((Y - X %*% B) %*% as.matrix(solve(s.ybar.y_r, tol = 1e-30)) %*% t(Y - X %*% B))))

    # BIC para las "no relevantes"
    if (diagonal) k <- m * (d + 1) + m
    else k <- m * (d + 1) + (m * (m + 1)) / 2
    bic.nc[b] <- 2 * loglik.nc - k * log(n)

    # BIC total
    Bic[b] <- bic.c[b] + bic.nc[b]
    if (!is.null(true.cl)) Ari[b] <- adjustedRandIndex(out$classification, true.cl)

    if (verb) print(b)
  }

  # Ensemble igual que antes
  print(Bic)
  cl.ens.1 <- data.frame(cl.m[, order(Bic, decreasing = TRUE)[1:B.star]])
  cl.ens.1.2 <- lapply(cl.ens.1, function(x) as.cl_membership(x))
  cl.consensus <- apply(cl_consensus(cl.ens.1.2, method = ensmethod)$.Data, 1, which.max)

  ari <- NULL
  if (!is.null(true.cl)) {
    ari <- adjustedRandIndex(cl.consensus, true.cl)
    names(ari) <- paste0("B.star=", B.star)
  }
  ensemble <- list(ari = ari, label.vec = cl.consensus)
  individual <- list(label.vec = cl.m, ari = Ari, bic = Bic, bic.GMM = bic.c, bic.reg = bic.nc)
  return(list(ensemble = ensemble, individual = individual))
}











