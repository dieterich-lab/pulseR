context("Parameters set up")

test_that("setting options work", {
  expect_error(fittingOptions(verbose = "x"))
  expect_error(fittingOptions(cores = -1))
  expect_equal(fittingOptions(options = list()), pulseR:::.defaultParams)
  
})

test_that("parameter initialization works", {
  expect_silent(setBoundaries(params = list(a = c(1, 2), b = c(1, 2)),
                              shared = list(d = c(1, 100))))
  # messed boundaries
  expect_error(setBoundaries(params = list(a = c(1, -2), b = c(1, 2)),
                             shared = list(d = c(1, 100))))
  # wrong length
  expect_error(setBoundaries(params = list(a = c(2), b = c(1, 2)),
                             shared = list(d = c(1, 100))))
  
  expect_silent(setTolerance(params = 1e-2))
  expect_error(setTolerance(params = -1e-2))
  expect_error(setTolerance(params = c(1, 2)))
  
  N <- 10
  pd <- list(count_data = matrix(ncol = 5, nrow = N))
  opts <- setBoundaries(params = list(a = c(1, 2), b = c(1, 4)),
                        shared = list(d = c(5, 15)))
  
  expect_silent(initParams(pd, opts))
  expect_silent(initParams(pd, opts, params = list(a = 1, b = 2)))
  # wrong name
  expect_error(initParams(pd, opts, shared = list(Y = 1)),
               regexp = "are not set")
  # forgot a param
  expect_error(initParams(pd, opts, params = list(a = rep(1,N))),
               regexp = "not equal")
  # out of range
  expect_error(initParams(pd, opts, params = list(a = -11, b = 2)),
               regexp = "within")
  # wrong number of genes
  expect_error(initParams(pd, opts, params = list(a = 1, b = 2)),
               regexp = "Length")
  
})
