## Tests for data with several time points
context("fitting with fraction factors for time dependent data")
set.seed(259)

nGenes <- 40
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
known <- addKnownToFormulas(formulas, formulaIndexes, conditions)

# create norm factors as 1...13, 13 is for the total ("A")
normFactors <- known$formulaIndexes[unique(names(known$formulaIndexes))]
normFactors <- normFactors[-grep("A", names(normFactors))]
normFactors <- c(list(total = 1), normFactors)
normFactors <- relist(seq_along(unlist(normFactors)), normFactors)

fractions <- as.character(interaction(conditions))
fractions[grep("A", fractions)] <- "total"

par <- list(size = 1e2)
par$alpha <-  .5
par <-  c(par, list(
  a = (1:nGenes) * 1e5, b = runif( nGenes,.1,.9)))
par$size <- 100000

allNormFactors <- multiplyList(normFactors, fractions)

counts <- generateTestDataFrom(
  formulas, formulaIndexes, allNormFactors, par, conditions)

pd <- PulseData(
  counts = counts,
  conditions = conditions,
  formulas = formulas,
  formulaIndexes = formulaIndexes,
  groups = fractions
)


#getNormIndex <- function(formulaIndexes, 
                

opts <- list()
opts$lb <- list(a=.1, b=.01)
opts$lb <- pulseR:::.b(opts$lb, par)
opts$ub <- list(a=1e7, b=.99)
opts$ub <- pulseR:::.b(opts$ub, par)
opts$lb$alpha <- .1
opts$ub$alpha <- 10
opts$lb$size <- 1
opts$ub$size <- 1e6

opts$lb$normFactors <- pulseR:::assignList(normFactors, .01)
opts$ub$normFactors <- pulseR:::assignList(normFactors, 20)
opts <- setTolerance(.1,shared = .1,fraction_factors = .1,options = opts)

opts$verbose <- "silent"
par$normFactors <- normFactors 

err <- function(x,y){
  vapply(intersect(names(x), names(y)), function(nx)
    max(1 - abs(x[[nx]])/abs(y[[nx]])), double(1))
}

test_that("gene params fitting works (together)", {
  skip("not testing gene params together")
  par2 <- par
  toOptimise <- c("a", "b")
  par2[toOptimise] <- opts$lb[toOptimise]
  res <- pulseR:::fitParams(pd, par, toOptimise, opts)
  expect_lt(max(err(res, par)), .1)
})

test_that("gene params fitting works (separately)", {
  par2 <- par
  toOptimise <- c("a", "b")
  par2[toOptimise] <- opts$lb[toOptimise]
  res <- pulseR:::fitParamsSeparately(pd, par, toOptimise, opts)
  expect_lt(max(err(res, par)), .1)
})

test_that("shared params fitting works", {
  par2 <- par
  toOptimise <- c("alpha")
  par2[toOptimise] <- opts$lb[toOptimise]
  res <- pulseR:::fitParams(pd, par, toOptimise, opts)
  expect_lt(max(err(res, par)), .1)
})

test_that("norm factors fitting works", {
  par2 <- par
  par2$normFactors <- pulseR:::assignList(par2$normFactors, 2)
  res <- pulseR:::fitNormFactors(pd, par2, opts)
  expect_lt(max(1-unlist(res)/unlist(par$normFactors)), .1)
})


test_that("all together fitting works", {
  par2 <- par
  par2["alpha"] <- .5
  par2[c("a", "b")] <- opts$ub[c("a", "b")]
  par2$normFactors <- pulseR:::assignList(par2$normFactors, 2)
  par2$a <- par$a
  opts$verbose <- "verbose"
  res <- pulseR:::fitModel(pd, par2, opts)
  expect_lt(max(1-unlist(res)/unlist(par)), .1)
})
  