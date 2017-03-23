
set.seed(259)

nGenes <- 10
nReplicates <- 3
nTime <- 3

formulas <- MeanFormulas(
  A = a,
  B =  a * b ^ time,
  C = alpha * a * (1 - b ^ time))

formulaIndexes <- list(
  A_fraction = 'A',
  B_fraction = c('B', 'C'),
  C_fraction = c('B', 'C'))


conditions <- data.frame(condition = rep(names(formulaIndexes), each = nTime),
                         time = rep(1:nTime, length(formulas) * nReplicates))
rownames(conditions) <- paste0("sample_", seq_along(conditions$condition))
known <- addKnownToFormulas(formulas, formulaIndexes, conditions)

# create norm factors as 1...13
normFactors <- known$formulaIndexes[unique(names(known$formulaIndexes))]
normFactors <- relist(seq_along(unlist(normFactors)), normFactors)

allNormFactors <- multiplyList(normFactors, names(known$formulaIndexes))
normFactors <- list(
  A_fraction = 1,
  B_fraction = c(2, .1),
  C_fraction = c(.1, 3))
allNormFactors <- multiplyList(normFactors, conditions[,1])

par <- list(size = 1e2)
par <-  c(par, list(
  a = (1:nGenes) * 1e5, b = runif( nGenes,.1,1)))
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
         rnbinom(n = length(res), mu = res, size = par$size)
       })
spikes <- do.call(cbind, spikes)
rownames(spikes) <- paste("spikes", seq_len(dim(spikes)[1]))
BSpikes <- rownames(spikes)[seq_len(numSpikes)]
CSpikes <- rownames(spikes)[numSpikes + seq_len(numSpikes)]

refGroup <- "A_fraction"

spikeLists <- list(
  A_fraction = list(c(BSpikes, CSpikes)),
  B_fraction = list(BSpikes, CSpikes),
  C_fraction = list(BSpikes, CSpikes))

counts <- rbind(counts, spikes)
spikeins <- list(refGroup = refGroup,
                 spikeLists = spikeLists)
pulseRSpikeinsData <- list(
  formulas = formulas,
  counts = counts,
  conditions = conditions,
  spikeins = spikeins,
  allNormFactors = allNormFactors,
  formulaIndexes = formulaIndexes,
  par = par
)
devtools::use_data(pulseRSpikeinsData, pkg ="pkg", overwrite = TRUE)