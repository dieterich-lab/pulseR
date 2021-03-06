## Tests for data with several time points
context("fitting with fraction factors for time dependent data")
set.seed(259)

nGenes <- 10
nReplicates <- 3
nTime <- 3

options <- setBoundaries(list(
  a = c(1, 1e10),
  b = c(.01, 1),
  alpha = c(.10, 100)
))

formulas <- MeanFormulas(
  A = a,
  B =  a * b ^ time,
  C = alpha * a * (1 - b ^ time))

formulaIndexes <- list(
  A_samp = 'A',
  B_samp = c('B', 'C'),
  C_samp = c('B', 'C'))


conditions <- data.frame(condition = rep(names(formulaIndexes), each = nTime),
                         time = rep(1:nTime, length(formulas) * nReplicates))
rownames(conditions) <- paste0("sample_", seq_along(conditions$condition))
known <- addKnownToFormulas(formulas, formulaIndexes, conditions)

# create norm factors as 1...13
normFactors <- known$formulaIndexes[unique(names(known$formulaIndexes))]
normFactors <- utils::relist(seq_along(unlist(normFactors)), normFactors)

allNormFactors <- multiplyList(normFactors, names(known$formulaIndexes))
normFactors <- list(
  A_samp = 1,
  B_samp = c(2, .1),
  C_samp = c(.1, 3))
allNormFactors <- multiplyList(normFactors, conditions[,1])

par <- list(size = 1e2)
par <-  c(par, list(a = runif(nGenes,1,1e5), b = runif( nGenes,.1,.9)))
par$alpha <- 5
par$size <- 100000


counts <- generateTestDataFrom(
  formulas, formulaIndexes, allNormFactors, par, conditions)

## make spikeins

numSpikes <- 10
spikeinsMeans <- runif(numSpikes, 10,10000)
spikes <- lapply(allNormFactors,
       function(x) {
         if (length(x) == 1)
           res <- rep(spikeinsMeans, 2) * x
         else
           res <- c(spikeinsMeans * x[1], spikeinsMeans * x[2])
         res
       })
spikes <- do.call(cbind, spikes)
rownames(spikes) <- paste("spikes", seq_len(dim(spikes)[1]))
Bpikes <- rownames(spikes)[seq_len(numSpikes)]
CSpikes <- rownames(spikes)[numSpikes + seq_len(numSpikes)]

refGroup <- "A_samp"

spikeLists <- list(
  A_samp = list(c(Bpikes, CSpikes)),
  B_samp = list(Bpikes, CSpikes),
  C_samp = list(Bpikes, CSpikes))

counts <- rbind(counts, spikes)
spikeins <- list(refGroup = refGroup,
                 spikeLists = spikeLists)
pd <- PulseData(
  counts = counts,
  conditions = conditions,
  formulas = formulas,
  formulaIndexes = formulaIndexes,
  spikeins = spikeins
)

options <- list()
options$lb <- list(a = .1, b = .01)
options$lb <- pulseR:::.b(options$lb, par)
options$ub <- list(a = 1e7, b = .99)
options$ub <- pulseR:::.b(options$ub, par)
options$lb$alpha <- .1
options$ub$alpha <- 100
options$lb$size <- 1
options$ub$size <- 1e6

options <- setTolerance(.01, shared = .01, options = options)

options$verbose <- "silent"

err <- function(x,y){
  vapply(intersect(names(x), names(y)), function(nx)
    max(1 - abs(x[[nx]])/abs(y[[nx]])), double(1))
}

test_that("gene params fitting works (together)", {
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
  res <- pulseR:::fitParamsSeparately(
    pd, par, toOptimise, knownGenePars = "p", options)
  expect_lt(max(err(res, par)), .1)
})

test_that("shared params fitting works", {
  par2 <- par
  toOptimise <- c("alpha")
  par2[toOptimise] <- options$lb[toOptimise]
  res <- pulseR:::fitParams(pd, par, toOptimise, options)
  expect_lt(max(err(res, par)), .1)
})
system.time(
test_that("all together fitting works", {
  par2 <- par
  for (p in c("a", "b")) {
    par2[[p]] <-
      runif(length(par[[p]]), options$lb[[p]], options$ub[[p]])
  }
  par2$alpha <- .1
  res <- pulseR:::fitModel(pd, par2, options)
  res$size <- NULL
  par$size <- NULL
  expect_lt(max(1 - unlist(res) / unlist(par)), .1)
})
) 
