#' Diverse Clustering Ensemble Plus
#'
#' Extension of the `dice()` function with future enhancements,
#' such as feature selection and dimensionality reduction.
#'
#' Currently, `dice_plus()` simply calls `dice()` with the same parameters.
#'
#' @inheritParams dice
#' @return A list with the following elements
#' \item{E}{raw clustering ensemble object}
#' \item{Eknn}{clustering ensemble object with knn imputation used on `E`}
#' \item{Ecomp}{flattened ensemble object with remaining missing entries imputed
#' by majority voting}
#' \item{clusters}{final clustering assignment from the diverse clustering
#' ensemble method}
#' \item{indices}{if `evaluate = TRUE`, shows cluster evaluation indices;
#' otherwise `NULL`}
#' @author Antonio C. Rodriguez
#' @export
#' @examples
#' data(hgsc)
#' dat <- hgsc[1:100, 1:50]
#' ref.cl <- strsplit(rownames(dat), "_") |>
#'   purrr::map_chr(2) |>
#'   factor() |>
#'   as.integer()
#' dice.obj <- dice_plus(dat, nk = 4, reps = 5, algorithms = "hc", cons.funs =
#' "kmodes", ref.cl = ref.cl, progress = FALSE, verbose = FALSE)
#' str(dice.obj, max.level = 2)
dice_plus <- function(...) {
dice(...)
}
