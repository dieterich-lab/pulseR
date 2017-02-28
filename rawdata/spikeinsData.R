
set.seed(259)

nGenes <- 10
nReplicates <- 3
nTime <- 3

options <- setBoundaries(
  params = list(a = c(1, 1e10), b = c(.01, 1)),
  shared = list(alpha = c(.10, 100)))

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
normFactors <- relist(seq_along(unlist(normFactors)), normFactors)

allNormFactors <- multiplyList(normFactors, names(known$formulaIndexes))
normFactors <- list(
  A_samp = 1,
  B_samp = c(2, .1),
  C_samp = c(.1, 3))
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