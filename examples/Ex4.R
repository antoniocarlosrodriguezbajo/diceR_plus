library(RPEClust)
library(future.apply)
library(mclust)
library(clue)

RPGMMClu_parallel <- function(x, true.cl=NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, modelNames = NULL, diagonal = FALSE, ensmethod="DWH", seed = 101, verb = FALSE){

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

system.time(out.clu_p <- RPGMMClu_parallel(Meat$x, Meat$y, g=5, B=100, B.star=10, verb=TRUE))
