## Tests for data with several time points
context("fitting with fraction factors for time dependent data")
source("test-utils.R")
set.seed(259)

nGenes <- 400
nReplicates <- 3
nTime <- 4

options <- list(
  lower_boundary = c(1,1e-3),
  upper_boundary = c(1e10, 1),
  lower_boundary_size = 1,
  upper_boundary_size = 1e9,
  lower_boundary_shared = .10,
  upper_boundary_shared = 100, 
  cores = 2
)
options$parscales <- c(1e5,1)
formulas <- MeanFormulas(A = a, B =  a * b^time, C= 1e7 * alpha^time)
conditions <- data.frame(condition = rep(names(formulas), each = nTime),
                         time = rep(1:nTime, length(formulas) * nReplicates))
rownames(conditions) <- paste0("sample_", seq_along(conditions$condition))
t <- pulseR:::addKnownShared(formulas, conditions)
formulas_known <- t$formulas
conditions_known <- data.frame(condition = t$conditions)

par <- list(size = 1e2)
par$shared_params <- list(alpha = 1)
par$individual_params <- data.frame(a = (1:nGenes) * 1e5, b = runif(nGenes,.1,.8))
rownames(par$individual_params) <- paste0("gene_", 1:nGenes)

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
  par2$individual_params$a <- guess
  par2$individual_params$b <- .3
  fit <- pulseR:::fitIndividualParameters(pd, par2, options)
  expect_lt(max(abs(1 - fit / par$individual_params)), .2)
})

test_that("shared params fitting works", {
  par2 <- par
  par2$shared_params <- list(alpha = 10)
  fit <- pulseR:::fitSharedParameters(pd, par2, options)
  expect_lt(max(abs(1 - unlist(fit) / unlist(par$shared_params))), .2)
})

test_that("overdispresion fitting works", {
  par2 <- par
  par2$size <- 1e4
  fit <- pulseR:::fitDispersion(pd, par2, options)
  expect_lt(max(abs(1 - unlist(fit) / unlist(par$size))), .2)
})
