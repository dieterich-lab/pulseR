## Tests for data with several time points
context("fraction normalisation for time dependent data")
source("test-utils.R")
set.seed(259)

nGenes <- 3
nReplicates <- 5
nTime <- 2

options <- list(
  lower_boundary = rep(1e-9, 2),
  upper_boundary = c(1e9, 1e9) - 1e-1,
  lower_boundary_size = 0,
  upper_boundary_size = 1e3,
  cores = 2
)
options$parscales <- mapply(max,
                            abs(options$upper_boundary),
                            abs(options$lower_boundary))

formulas <- MeanFormulas(A = a, B =  b * time)
conditions <- data.frame(condition = rep(c("A","B"), each=nTime),
                         time = rep(1:nTime, 2 * nReplicates))
rownames(conditions) <- paste0("sample_", seq_along(conditions$condition))
t <- pulseR:::addKnownShared(formulas, conditions)
formulas_known <- t$formulas
conditions_known <- data.frame(condition=t$conditions)

par <- list(size=1e7)
fractions <- factor(conditions_known$condition)
par$individual_params <- data.frame(a=1:nGenes * 1e7, b=rep(5e7, nGenes))
rownames(par$individual_params) <- paste0("gene_", 1:nGenes)

par$fraction_factors <- 1 * (1:(length(levels(fractions))-1))
#par$fraction_factors <- rep(2,length(levels(fractions))-1)
#par$fraction_factors <- rep(1:5,2)[-1]
counts <- generateTestDataFrom(formulas_known,par,conditions_known, fractions)

pd <- PulseData(
    count_data = counts,
    conditions = conditions,
    formulas   = formulas,
    fractions  = ~condition+time)
normalise(pd) 
par2 <- par
par2$individual_params$a <- 1e3
par2$individual_params$b <- 1e3
test_that("individual parameters fitting works", {
  fit <- pulseR:::fitIndividualParameters(pd, par, options)
  expect_lt(max(abs(1 - fit / par$individual_params)), .2)
})
#fit <- fitModel(pd,par2,options)

#test_that("fitting works for time-series", {
#  expect_gt(.3,
#                   max(abs((fit$par$individual_params - par$individual_params) /
#                             par$individual_params
#                   )))
#})