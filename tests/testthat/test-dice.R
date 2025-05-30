data(hgsc)
ref.cl <- rownames(hgsc) |>
  strsplit("_") |>
  vapply(`[`, 2, FUN.VALUE = character(1)) |>
  factor() |>
  as.integer()

test_that("dice works with one algorithm, one consensus funs", {
  dice.obj <- dice(
    hgsc,
    nk = 4,
    algorithms = "hc",
    cons.funs = "kmodes",
    progress = FALSE,
    verbose = FALSE
  )
  expect_length(dice.obj, 5)
  expect_equal(dim(dice.obj$clusters), c(nrow(hgsc), 1))
})

test_that("dice works with multiple algorithms, consensus funs, trimming, and
          reference class", {
            dice.obj <- dice(
              hgsc,
              nk = 4,
              reps = 5,
              algorithms = c("hc", "diana"),
              cons.funs = c("kmodes", "majority"),
              trim = TRUE,
              n = 2,
              ref.cl = ref.cl,
              progress = FALSE,
              verbose = FALSE
            )
  expect_length(dice.obj, 5)
  expect_true(is.matrix(dice.obj$clusters))
})

test_that("single algorithm and single consensus return same results", {
  dice.obj <- dice(
    hgsc,
    nk = 4,
    reps = 5,
    algorithms = "km",
    cons.funs = "CSPA",
    ref.cl = ref.cl,
    progress = FALSE,
    verbose = FALSE
  )
  ind.obj <- dice.obj$indices$ei$`4`
  expect_equal(unname(unlist(ind.obj[1, -1])), unname(unlist(ind.obj[2, -1])))
})

test_that("indices slot returns NULL if evaluate specified as FALSE", {
  dice.obj <- dice(
    hgsc,
    nk = 4,
    reps = 3,
    algorithms = "hc",
    cons.funs = "kmodes",
    ref.cl = ref.cl,
    evaluate = FALSE,
    progress = FALSE,
    verbose = FALSE
  )
  expect_null(dice.obj$indices)
})

test_that("relabelling uses 1st col if more than 1 cons.funs and no ref.cl", {
  dice.obj <- dice(
    hgsc,
    nk = 4,
    reps = 3,
    algorithms = "hc",
    cons.funs = c("kmodes", "majority"),
    evaluate = FALSE,
    progress = FALSE,
    verbose = FALSE
  )
  expect_error(dice.obj, NA)
})

test_that("cluster size prepended when multiple k requested", {
  dice.obj <- dice(
    hgsc,
    nk = 3:4,
    reps = 3,
    algorithms = "hc",
    cons.funs = "kmodes",
    k.method = "all",
    evaluate = FALSE,
    progress = FALSE,
    verbose = FALSE
  )
  expect_true(all(grepl("k=", colnames(dice.obj$clusters))))
})

test_that("algorithm vs internal index heatmap works", {
  dice.obj <- dice(
    hgsc,
    nk = 4,
    reps = 3,
    algorithms = "hc",
    cons.funs = "kmodes",
    ref.cl = ref.cl,
    evaluate = FALSE,
    plot = TRUE,
    progress = FALSE,
    verbose = FALSE
  )
  expect_error(dice.obj, NA)
})
