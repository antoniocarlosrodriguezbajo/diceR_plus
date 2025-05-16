#'
#' @title Generation of random matrices
#'
#' @description This function generates \code{B} random matrices of dimension \code{p} by \code{d} by using the Haar measure.
#'
#' @param p The original number of variables.
#' @param d The reduced dimension.
#' @param B The number of projections.
#' @return A single matrix of dimension \code{p} by \code{d*B} containing \code{B} random matrices of dimension \code{p} by \code{d}.
#' @references Anderlucci, Fortunato, Montanari (2019) <arXiv:1909.10832>
#' @examples
#' R<-generateRP(p=100,d=2,B=10)
#' dim(R)
#' @export



generateRP<-function (p, d, B){
  if (p < d) stop("d must be less than p")
  R0 <- matrix(1/sqrt(p) * rnorm(p * d * B, 0, 1), p, d * B)
  Q <- matrix(sapply(1:B - 1, FUN = function(s) {
    qr.Q(qr(R0[, s * d + 1:d]))[, 1:d]}), p, B * d)
  return(Q)}
