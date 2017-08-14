
set.seed(259)

nGenes <- 50
nReplicates <- 2
nTime <- 3

formulas <- MeanFormulas(
  total = mu,
  labelled =  mu * exp(- d* time),
  unlabelled =  mu *exp(-d*time))

formulaIndexes <- list(total_fraction = "total",
                       pull_down = c("labelled", "unlabelled"))

conditions <- data.frame(
  fraction = rep(names(formulaIndexes), each = nTime),
  time = rep(1:nTime, length(formulaIndexes) * nReplicates)
)
rownames(conditions) <- paste0("sample_", seq_along(conditions$fraction))
known <- addKnownToFormulas(formulas, formulaIndexes, conditions)

# create norm factors as 1...13
normFactors <- known$formulaIndexes[unique(names(known$formulaIndexes))]
normFactors <- relist(seq_along(unlist(normFactors)), normFactors)

allNormFactors <- multiplyList(normFactors, names(known$formulaIndexes))
normFactors <- list(
  total_fraction = 1,
  pull_down = c(2, .1))
allNormFactors <- multiplyList(normFactors, conditions[,1])

par <- list(size = 1e2)
par <-  c(par, list(
  mu = (1:nGenes) * 1e4, d = runif( nGenes,.1,.3)))
par$size <- 10000


counts <- generateTestDataFrom(
  formulas, formulaIndexes, allNormFactors, par, conditions)

## make spikeins

numSpikes <- 5
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
labelled <- rownames(spikes)[seq_len(numSpikes)]
unlabelled <- rownames(spikes)[numSpikes + seq_len(numSpikes)]

refGroup <- "total_fraction"

spikeLists <- list(
  total_fraction = list(c(labelled, unlabelled)),
  pull_down = list(labelled, unlabelled))

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
