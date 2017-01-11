## Tests for data with several time points
context("time dependent data")
source("test-utils.R")
set.seed(259)

options <- list(
  lower_boundary = rep(1e-9, 2),
  upper_boundary = c(1e9, 1e9) - 1e-1,
  lower_boundary_size = 0,
  upper_boundary_size = 1e3,
  lower_boundary_fraction = .1,
  upper_boundary_fraction = 10,
  cores = 2
)
options$parscales <- mapply(max,
                            abs(options$upper_boundary),
                            abs(options$lower_boundary))

formulas <- MeanFormulas(A = a, B = a + b * time)
conditions <- data.frame(condition = c(rep("A", 5), rep("B", 5)),
                         time = rep(1:5, 2))
rownames(conditions) <-
  paste0("sample_", seq_along(conditions$condition))

par <- list(size=100)
par$individual_params <- data.frame(a=c(1e6,1e7), b=c(1,3)*1e6)
rownames(par$individual_params) <-
  paste0("gene_", 1:length(par$individual_params))
counts <- generateTestDataFrom(formulas,par,conditions)

pd <- PulseData(
    count_data = counts,
    conditions = conditions,
    formulas   = formulas,
    fractions  = ~condition+time)
normalise(pd) 
par2 <- par
par2$individual_params$a <- 1e3
par2$individual_params$b <- 1e3
par2$fraction_factors <- rep(1, length(pd$fraction)-1)
fit <- fitModel(pd,par2,options)
test_that("fitting works for time-series", {
  expect_gt(.3,
                   max(abs((fit$par$individual_params - par$individual_params) /
                             par$individual_params
                   )))
})