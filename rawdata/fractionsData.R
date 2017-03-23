
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

pulseRFractionData <- list(
  formulas = formulas,
  counts = counts,
  conditions = conditions,
  fractions = fractions,
  formulaIndexes = formulaIndexes,
  normFactors = normFactors,
  par = par
)
devtools::use_data(pulseRFractionData, overwrite = TRUE)
