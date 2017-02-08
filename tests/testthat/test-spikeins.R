## Tests for data with several time points
context("fitting with fraction factors for time dependent data")
source("test-utils.R")
set.seed(259)

nGenes <- 400
nReplicates <- 3
nTime <- 4

options <- setBoundaries(params = list(lb = c(1, 1e-3),
                                       up = c(1e10, 1)),
                         shared = list(lb = .10,
                                       ub = 100))

options$cores <- 1

formulas <- MeanFormulas(A = a * p, B =  a * b^time, C= 1e7 * alpha^time)
conditions <- data.frame(condition = rep(names(formulas), each = nTime),
                         time = rep(1:nTime, length(formulas) * nReplicates))
rownames(conditions) <- paste0("sample_", seq_along(conditions$condition))
t <- pulseR:::addKnownShared(formulas, conditions)
formulas_known <- t$formulas
conditions_known <- data.frame(condition = t$conditions)

par <- list(size = 1e2)
par$known <-data.frame(p=1:nGenes + 1)
par$shared <- list(alpha = 1)
par$params <- data.frame(a = (1:nGenes) * 1e5, b = runif(nGenes,.1,.8))
rownames(par$params) <- paste0("gene_", 1:nGenes)

counts <- generateTestDataFrom(formulas_known,par,conditions_known)
# make mock spikeins
nSpikeIns <- 10
spikeins <- rnbinom(mu=1e5, size = par$size,n = ncol(counts) * nSpikeIns) 
spikeins <- matrix(spikeins, ncol=ncol(counts))
rownames(spikeins) <- paste("spike", 1:nSpikeIns, sep="_")
counts <- rbind(counts, spikeins)
pd <- PulseData(
    count_data = counts,
    conditions = conditions,
    formulas   = formulas,
    spikeins = rownames(spikeins))

test_that("individual parameters fitting works", {
  par2 <- par
  guess <-  apply(pd$count_data[, conditions$condition == "A"], 1, mean)
  par2$params$a <- guess
  par2$params$b <- .3
  fit <- pulseR:::fitGeneParameters(pd, par2, options)
  expect_lt(max(abs(1 - fit / par$params)), .2)
})

test_that("shared params fitting works", {
  par2 <- par
  par2$shared_params <- list(alpha = 10)
  fit <- pulseR:::fitSharedParameters(pd, par2, options)
  expect_lt(max(abs(1 - unlist(fit) / unlist(par$shared))), .2)
})

test_that("overdispresion fitting works", {
  par2 <- par
  par2$size <- 1e4
  fit <- pulseR:::fitDispersion(pd, par2, options)
  expect_lt(max(abs(1 - unlist(fit) / unlist(par$size))), .2)
})

