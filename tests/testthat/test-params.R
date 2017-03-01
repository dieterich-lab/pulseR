context("Parameters set up")

test_that("setting options work", {
  expect_error(setFittingOptions(verbose = "x"))
  expect_error(setFittingOptions(cores = -1))
  expect_equal(setFittingOptions(options = list()), pulseR:::.defaultParams)
  
})

test_that("parameter initialization works", {
  skip("TODO")
  expect_silent(setBoundaries(list(a = c(1, 2), b = c(1, 2), d = c(1, 100))))
  # messed boundaries
  expect_error(setBoundaries(list(a = c(1, -2), b = c(1, 2), d = c(1, 100))))
  # wrong length
  expect_error(setBoundaries(list(a = c(2), b = c(1, 2), d = c(1, 100))))
  
  expect_silent(setTolerance(params = 1e-2))
  expect_error(setTolerance(params = -1e-2))
  expect_error(setTolerance(params = c(1, 2)))
  
})
