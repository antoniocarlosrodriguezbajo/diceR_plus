#' @title Generation of random projection matrices
#'
#' @description This function generates \code{B} random projection matrices of dimension \code{p} by \code{d}, either using orthonormal columns from the Haar measure (via Gaussian QR decomposition), or using sparse Achlioptas-type projections.
#'
#' @param p The original number of variables (input dimension).
#' @param d The reduced target dimension.
#' @param B The number of projections to generate (default is 1).
#' @param method_rp Character string indicating the type of random projection: \code{"gaussian"} for Haar-orthonormal projection matrices, or \code{"achlioptas"} for sparse projections (default is \code{"gaussian"}).
#' @param s Sparsity parameter used in the \code{"achlioptas"} method. Larger values yield sparser matrices (default is 3).
#'
#' @return A numeric matrix of size \code{p} by \code{d * B}, containing \code{B} random projection matrices concatenated column-wise.
#'
#' @references
#' Anderlucci, L., Fortunato, L., & Montanari, A. (2019). A general framework for dimension reduction in multivariate and functional data. \emph{arXiv preprint} arXiv:1909.10832.
#' Achlioptas, D. (2003). Database-friendly random projections: Johnson–Lindenstrauss with binary coins. \emph{J. Comput. System Sci.}, 66(4), 671–687.
#'
#' @examples
#' # Generate 10 Gaussian orthonormal projection matrices
#' G <- generateRP(p = 100, d = 5, B = 10, method_rp = "gaussian")
#' dim(G)
#'
#' # Generate 5 sparse Achlioptas-type projection matrices with default sparsity
#' A <- generateRP(p = 100, d = 5, B = 5, method_rp = "achlioptas", s = 3)
#' dim(A)
#'
#' @export

generateRP <- function(p, d, B = 1, method_rp = c("gaussian", "achlioptas"), s = 3) {
  method_rp <- match.arg(method_rp)
  if (p < d) stop("d must be less than or equal to p")

  if (method_rp == "gaussian") {
    R0 <- matrix(rnorm(p * d * B, mean = 0, sd = 1) / sqrt(p), nrow = p)
    Q_blocks <- lapply(0:(B - 1), function(s_idx) {
      cols <- (s_idx * d + 1):(s_idx * d + d)
      qr.Q(qr(R0[, cols]))[, 1:d]
    })
    Q <- do.call(cbind, Q_blocks)

  } else if (method_rp == "achlioptas") {
    prob <- 1 / (2 * s)
    values <- c(sqrt(s), 0, -sqrt(s))
    probs <- c(prob, 1 - 2 * prob, prob)

    R_blocks <- lapply(1:B, function(b) {
      R <- matrix(sample(values, size = p * d, replace = TRUE, prob = probs), nrow = p)
      R / sqrt(p)
    })
    Q <- do.call(cbind, R_blocks)
  }

  return(Q)
}

