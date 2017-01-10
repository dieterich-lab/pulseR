## Tests for data with several time points
context("fraction normalisation for time dependent data")
source("test-utils.R")
set.seed(259)

nGenes <- 20
nReplicates <- 3
nTime <- 7

options <- list(
  lower_boundary = c(1,1e-3),
  upper_boundary = c(1e10, 1),
  lower_boundary_size = 1,
  upper_boundary_size = 1e9,
  lower_boundary_shared = .10,
  upper_boundary_shared = 100, 
  cores = 2
)
options$parscales <- mapply(max,
                            abs(options$upper_boundary),
                            abs(options$lower_boundary))

formulas <- MeanFormulas(A = a, B =  a * b^time, C= alpha * a)
conditions <- data.frame(condition = rep(c("A", "B", "C"), each = nTime),
                         time = rep(1:nTime, length(formulas) * nReplicates))
rownames(conditions) <- paste0("sample_", seq_along(conditions$condition))
t <- pulseR:::addKnownShared(formulas, conditions)
formulas_known <- t$formulas
conditions_known <- data.frame(condition = t$conditions)

par <- list(size = 1e7)
par$shared_params <- list(alpha = 1)
fractions <- factor(conditions_known$condition)
par$individual_params <-
  data.frame(a = (1:nGenes) * 1e7, b = rep(.8, nGenes))
rownames(par$individual_params) <- paste0("gene_", 1:nGenes)

par$fraction_factors <- 1 * (1:(length(levels(fractions)) - 1))
#par$fraction_factors <- rep(1,length(levels(fractions))-1)
#par$fraction_factors <- rep(1:5,2)[-1]
counts <- generateTestDataFrom(formulas_known,par,conditions_known, fractions)

pd <- PulseData(
    count_data = counts,
    conditions = conditions,
    formulas   = formulas,
    fractions  = ~condition+time)
normalise(pd) 
par2 <- par
par2$individual_params$a <- 5e7
par2$individual_params$b <- .2
test_that("individual parameters fitting works", {
  fit <- pulseR:::fitIndividualParameters(pd, par2, options)
  expect_lt(max(abs(1 - fit / par$individual_params)), .2)
})

test_that("shared params fitting works", {
  options$lower_boundary_shared <- .10
  options$upper_boundary_shared <- 100
  par2 <- par
  par2$shared_params <- list(alpha = 10)
  fit <- pulseR:::fitSharedParameters(pd, par2, options)
  expect_lt(max(abs(1 - unlist(fit) / unlist(par$shared_params))), .2)
})

test_that("overexpression fitting works", {
  par2 <- par
  par2$size <- 1e4
  fit <- pulseR:::fitDispersion(pd, par2, options)
  expect_lt(max(abs(1 - unlist(fit) / unlist(par$size))), .2)
})

test_that("fraction factors fitting works", {
  par2 <- par
  par2$fraction_factors <- rep(10, length(par$fraction_factors))
  options$lower_boundary_fraction <- rep(.1, length(par$fraction_factors))
  options$upper_boundary_fraction <- rep(100, length(par$fraction_factors))
  fit <- pulseR:::fitFractions(pd, par2, options)
  expect_lt(max(abs(1 - unlist(fit) / unlist(par$fraction_factors))), .2)
})

test_that("all together fitting works", {
  par2 <- par
  par2$size <- 1e4
  
  options$lower_boundary_shared <- .10
  options$upper_boundary_shared <- 100
  
  par2$fraction_factors <- rep(10, length(par$fraction_factors))
  options$lower_boundary_fraction <- rep(.1, length(par$fraction_factors))
  options$upper_boundary_fraction <- rep(10, length(par$fraction_factors))
  #options$verbose <- "verbose"
  fit <- fitModel(pd, par2, options)
  expect_gt(.3,
                   max(abs((fit$par$individual_params - par$individual_params) /
                             par$individual_params
                   )))
})
