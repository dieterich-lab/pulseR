## Tests for data with several time points
context("fitting with fraction factors for time dependent data")
source("test-utils.R")
set.seed(259)

nGenes <- 40
nReplicates <- 3
nTime <- 4

options <- setBoundaries(
  params = list(a = c(1, 1e10), b = c(.01, 1)),
  shared = list(alpha = c(.10, 100)),
  fraction_factors = list( c(.1,10)))

options$cores <- 1

formulas <- MeanFormulas(A = a, B =  a * b^time, C= 1e7 * alpha^time)
conditions <- data.frame(condition = rep(names(formulas), each = nTime),
                         time = rep(1:nTime, length(formulas) * nReplicates))
rownames(conditions) <- paste0("sample_", seq_along(conditions$condition))
t <- pulseR:::addKnownShared(formulas, conditions)
formulas_known <- t$formulas
conditions_known <- data.frame(condition = t$conditions)

par <- list(size = 1e2)
par$shared <- list(alpha = 1)
fractions <- factor(conditions_known$condition)
par$params<- data.frame(a = (1:nGenes) * 1e5, b = runif(nGenes,.1,.8))
rownames(par$params) <- paste0("gene_", 1:nGenes)

#par$fraction_factors <- 1 * (1:(length(levels(fractions)) - 1))
#par$fraction_factors <- rep(1,length(levels(fractions))-1)
par$fraction_factors <- runif(length(levels(fractions)), 1, 5)

par$fraction_factors[1] <- 1

counts <- generateTestDataFrom(formulas_known,par,conditions_known, fractions)

pd <- PulseData(
    count_data = counts,
    conditions = conditions,
    formulas   = formulas,
    fractions  = ~condition+time)
normalise(pd) 

test_that("gene-specific parameters fitting works", {
  par2 <- par
  guess <-  apply(counts[, conditions$condition == "A"], 1, mean)
  par2$params$a <- guess
  par2$params$b <- .3
  fit <- pulseR:::fitGeneParameters(pd, par2, options)
  expect_lt(max(abs(1 - fit / par$params)), .2)
})

test_that("shared params fitting works", {
  par2 <- par
  par2$shared <- list(alpha = 10)
  fit <- pulseR:::fitSharedParameters(pd, par2, options)
  expect_lt(max(abs(1 - unlist(fit) / unlist(par$shared))), .2)
})

test_that("overdispersion fitting works", {
  par2 <- par
  par2$size <- 1e4
  fit <- pulseR:::fitDispersion(pd, par2, options)
  expect_lt(max(abs(1 - unlist(fit) / unlist(par$size))), .2)
})

test_that("fraction factors fitting works", {
  par2 <- par
  par2$fraction_factors <- rep(1, length(par$fraction_factors))
  fit <- pulseR:::fitFractions(pd, par2, options)
  expect_lt(max(abs(1 - unlist(fit) / unlist(par$fraction_factors))), .2)
})
