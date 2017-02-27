context("PulseData creation and normalisation")

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

test_that( "pd can be created without groups", {
    pd <- PulseData(
      counts = counts,
      conditions = conditions,
      formulas = formulas,
      formulaIndexes = formulaIndexes
    )
})

test_that( "groups can be set by vector", {
  expect_silent({
    pd <- PulseData(
      counts = counts,
      conditions = conditions,
      formulas = formulas,
      formulaIndexes = formulaIndexes,
      groups = conditions[,1]
    )
  })
})

test_that( "groups can be set by formula", {
  expect_silent({
    pd <- PulseData(
      counts = counts,
      conditions = conditions,
      formulas = formulas,
      formulaIndexes = formulaIndexes,
      groups = ~ time
    )
  })
})

