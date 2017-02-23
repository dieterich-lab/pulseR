## Tests for data with several time points
context("fitting with fraction factors for time dependent data")
set.seed(259)

nGenes <- 4
nReplicates <- 3
nTime <- 4


formulas <- MeanFormulas(A = a, B =  a * b^time, C = a * alpha^time)


formulaIndexes <- list(
  A_samp = 'A',
  B_samp = c('B', 'A'),
  C_samp = c('C'))

normFactors <- list(
  A_samp = c(1),
  B_samp = c(1, .1),
  C_samp = 2
)

conditions <- data.frame(condition = rep(names(formulaIndexes), each = nTime),
                         time = rep(1:nTime, length(formulas) * nReplicates))
rownames(conditions) <- paste0("sample_", seq_along(conditions$condition))

normFactors <- pulseR:::multiplyList(normFactors,conditions[,1])

par <- list(size = 1e2)
par$alpha <-  .5
par <-  c(par, list(
  a = (1:nGenes) * 1e5, b = runif( nGenes,.1,.9)))
par$size <- 100000



counts <- generateTestDataFrom(
  formulas, formulaIndexes, normFactors, par, conditions)

pd <- PulseData(
  counts = counts,
  conditions = conditions,
  formulas = formulas,
  formulaIndexes = formulaIndexes
)
                
opts <- list()
opts$lb <- list(a=.1, b=.01)
opts$lb <- pulseR:::.b(opts$lb, par)
opts$ub <- list(a=1e7, b=.99)
opts$ub <- pulseR:::.b(opts$ub, par)
opts$lb$alpha <- .1
opts$ub$alpha <- 10

par$normFactors <- normFactors
err <- function(x,y){
  vapply(intersect(names(x), names(y)), function(nx)
    max(1 - abs(x[[nx]])/abs(y[[nx]])), double(1))
}

test_that("gene params fitting works (together)", {
  par2 <- par
  toOptimise <- c("a", "b")
  par2[toOptimise] <- opts$lb[toOptimise]
  res <- fitParams(pd, par, toOptimise, opts)
  expect_lt(max(err(res, par)), .1)
})

test_that("gene params fitting works (separately)", {
  par2 <- par
  toOptimise <- c("a", "b")
  par2[toOptimise] <- opts$lb[toOptimise]
  res <- fitParamsSeparately(pd, par, toOptimise, opts)
  expect_lt(max(err(res, par)), .1)
})

test_that("shared params fitting works", {
  par2 <- par
  toOptimise <- c("alpha")
  par2[toOptimise] <- opts$lb[toOptimise]
  res <- fitParams(pd, par, toOptimise, opts)
  expect_lt(max(err(res, par)), .1)
})


