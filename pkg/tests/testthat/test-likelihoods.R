context("likelihood functions and formula manipulations")

test_that("correct fraction contamination",{
  formulas <- MeanFormulas(label = a,
                           unlabel = b,
                           extra = d + 1)
  expected <- quote((1 - p) * (a) + p * (b))
  expect_equal(contaminate(formulas, "label", "unlabel", "p"), expected)
})