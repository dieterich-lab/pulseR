## Tests for data with several time points
context("fitting with fraction factors for time dependent data")
set.seed(259)

nGenes <- 10
nReplicates <- 3
nTime <- 3


formulas <- MeanFormulas(A = a * p, B =  a * b^time, C = a * (1 - b^time))


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
normFactors[grep("B", names(normFactors))] <- list(c(3,.2))

fractions <- as.character(interaction(conditions))
fractions[grep("A", fractions)] <- "total"

par <- list(size = 1e2)
par <-  c(par, list(
  a = runif(nGenes, 1, 1e5), b = runif( nGenes,.1,.91)))
par$p <- runif(nGenes, 1, 2)
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


options <- list()
options$lb <- list(a = .1, b = .01)
options$lb <- pulseR:::.b(options$lb, par)
options$ub <- list(a = 1e7, b = .99)
options$ub <- pulseR:::.b(options$ub, par)
options$lb$size <- 1
options$ub$size <- 1e6

options$lb$normFactors <- pulseR:::assignList(normFactors, .01)
options$ub$normFactors <- pulseR:::assignList(normFactors, 20)
options <- setTolerance(.01,shared = .01, normFactors = .01,options = options)

options$verbose <- "silent"
par$normFactors <- normFactors 

err <- function(x,y){
  vapply(intersect(names(x), names(y)), function(nx)
    max(1 - abs(x[[nx]])/abs(y[[nx]])), double(1))
}

test_that("gene params fitting works (together)", {
  skip("skip long test")
  par2 <- par
  toOptimise <- c("a", "b")
  par2[toOptimise] <- options$lb[toOptimise]
  res <- pulseR:::fitParams(pd, par, toOptimise, options)
  expect_lt(max(err(res, par)), .1)
})

test_that("gene params fitting works (separately)", {
  par2 <- par
  toOptimise <- c("a", "b")
  par2[toOptimise] <- options$lb[toOptimise]
  res <- pulseR:::fitParamsSeparately(pd = pd,
                                 par = par,
                                 knownNames = "p",
                                 namesToOptimise = c("a", "b"), 
                                 options = options)
  expect_lt(max(err(res, par)), .1)
})

test_that("norm factors fitting works", {
  par2 <- par
  par2$normFactors <- pulseR:::assignList(par2$normFactors, 2)
  res <- pulseR:::fitNormFactors(pd, par2, options)
  expect_lt(max(1 - unlist(res) / unlist(par$normFactors)), .1)
})


test_that("all together fitting works", {
  par2 <- par
  for (p in c("a", "b")) {
    par2[[p]] <-
      runif(length(par[[p]]), options$lb[[p]], options$ub[[p]])
  }
  options$verbose <- "verbose"
  res <- pulseR:::fitModel(pd, par2, options)
  res$size <- NULL
  par$size <- NULL
  expect_lt(max(1 - unlist(res) / unlist(par)), .1)
})
  