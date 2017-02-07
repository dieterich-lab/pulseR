## Tests for data with several time points
context("fitModel with fractions")
source("test-utils.R")
set.seed(259)
nGenes <- 40
nReplicates <- 3
nTime <- 4

options <- list(
  lower_boundary = c(1,1e-3),
  upper_boundary = c(1e10, 1-1e-3),
  lower_boundary_size = 1,
  upper_boundary_size = 1e9,
  cores = 2
)
options$parscales <- c(1e5,1)

formulas <- MeanFormulas(A = a, B =  a * b^time, C= a * (1-b^time))
conditions <- data.frame(condition = rep(names(formulas), each = nTime),
                         time = rep(1:nTime, length(formulas) * nReplicates))
rownames(conditions) <- paste0("sample_", seq_along(conditions$condition))
t <- pulseR:::addKnownShared(formulas, conditions)
formulas_known <- t$formulas
conditions_known <- data.frame(condition = t$conditions)

par <- list(size = 1e2)
fractions <- factor(conditions_known$condition)
par$individual_params <-
  data.frame(a = (1:nGenes) * 1e5, b = runif( nGenes,.1,.8))
rownames(par$individual_params) <- paste0("gene_", 1:nGenes)

par$fraction_factors <- 1 * (1:(length(levels(fractions)) ))
#par$fraction_factors <- rep(1,length(levels(fractions))-1)
#par$fraction_factors <- runif(length(levels(fractions)) - 1, 1, 5)
options$lower_boundary_fraction <- rep(.1, length(par$fraction_factors))
options$upper_boundary_fraction <- rep(100, length(par$fraction_factors))

counts <- generateTestDataFrom(formulas_known,par,conditions_known, fractions)

pd <- PulseData(
    count_data = counts,
    conditions = conditions,
    formulas   = formulas,
    fractions  = ~condition+time)
normalise(pd) 

test_that("all together fitting works", {
  #skip("...")
  par2 <- par
  guess <-  apply(counts[, conditions$condition == "A"], 1, mean)
  par2$individual_params$a <- guess
  par2$individual_params$b <- runif(length(par$individual_params$a),.1,.8)
  par2$size <- 1e4
  par2$fraction_factors <- rep(1, length(par$fraction_factors))
  options$verbose <- "verbose"
  fit <- fitModel(pd, par2, options)
  expect_gt(.3,
                   max(abs((fit$par$individual_params - par$individual_params) /
                             par$individual_params
                   )))
  expect_lt(max(abs(
    1 - unlist(fit$par$fraction_factors) / unlist(par$fraction_factors)
  )), .3)
})
