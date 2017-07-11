context("Confidence intervals")
set.seed(259)



formulas <- MeanFormulas(X = mu, Y = nu)


formulaIndexes <- list(
  EX = 'X',
  EXandY = c('X', 'Y'))

normFactors <- list(
  EX = c(1),
  EXandY = c(1, .1)
)

nTime <- 1
nReplicates <- 4
conditions <- data.frame(condition = rep(names(formulaIndexes), each = nTime),
                         time = rep(1:nTime, length(formulas) * nReplicates))
rownames(conditions) <- paste0("sample_", seq_along(conditions$condition))
known <- addKnownToFormulas(formulas, formulaIndexes, conditions)
normFactors <- known$formulaIndexes[unique(names(known$formulaIndexes))]

fractions <- as.character(interaction(conditions))

nGenes <- 2
par <- list(size = 1e4)
par <-  c(par, list(
  mu = runif(nGenes, 100, 1000), nu = runif(nGenes,100,1000)))

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
options$lb <- list(mu = 1, nu = 1)
options$lb <- pulseR:::.b(options$lb, par)
options$ub <- list(mu = 1e4, nu = 1e4)
options$ub <- pulseR:::.b(options$ub, par)
options$lb$size <- 1
options$ub$size <- 1e6

options$lb$normFactors <- pulseR:::assignList(normFactors, .01)
options$ub$normFactors <- pulseR:::assignList(normFactors, 20)
options <- setTolerance(.01,shared = .01, normFactors = .01,options = options)

options$verbose <- "silent"
par$normFactors <- normFactors 
fit <- fitModel(pd, par,options)

test_that("plGene is zero at optimum", {
   expect_lte(
     abs(plGene("mu",1,fit, pd,options)(fit$mu[1])), 1e-4)
   expect_lte(
     abs(pl(list("mu",1),fit, pd,options)(fit$mu[1])), 1e-4)
})

test_that("profile estimations on the interval", {
   prof <- profile(list("mu", 1), pd, fit, options,
                   interval = rep(fit$mu[1], 2), numPoints = 1) 
   expect_lte(abs(prof$logL), 1e-6)
   prof <- profileOnlyGene("mu", 1, pd, fit, options,
                   interval = rep(fit$mu[1], 2), numPoints = 1) 
   expect_lte(abs(prof$logL), 1e-6)
})


test_that("ci calculation", {
  cis <- ciGene("mu",1,pd,fit,options)
  optimum <- evaluateLikelihood(fit, pd)
  vapply(cis, function(x) {
    p <- .assignElement(fit, list("mu",1), x)
    evaluateLikelihood(p, pd) - optimum
  }, double(1))
  cis <- ci(list("mu", 1), pd, fit, options) 
})