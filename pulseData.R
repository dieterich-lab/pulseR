
PulseData <- function(count_data,
                      conditions,
                      formulas,
                      spikeins = rownames(count_data)) {
  e <- new.env()
  e$count_data <- count_data
  e$conditions <- conditions
  e$formulas <- formulas
  e$spikeins <- spikeins
  e
}

findDeseqFactorsSingle <- function(count_data)
{
  loggeomeans <- rowMeans(log(count_data))
  deseqFactors <-  apply(count_data, 2, function(x) {
    exp(median(log(x) - loggeomeans, na.rm = TRUE))
  })
  deseqFactors
}

# conditions  - factor to split samples for normalisation
findDeseqFactors <- function(count_data,
                             conditions,
                             spikeins) {
  deseqFactors <- lapply(split(rownames(conditions), conditions),
                         function(samples) {
                           findDeseqFactorsSingle(count_data[spikeins, samples])
                         })
  names(deseqFactors) <- NULL
  unlist(deseqFactors)[colnames(count_data)]
}

normalise <- function(pulseData, fractions) {
  if (missing(fractions)) {
    conditions <- pulseData$conditions[, 1]
  } else {
    conditions <- pulseData$conditions[, all.vars(fractions)]
  }
  pulseData$norm_factors <- findDeseqFactors(pulseData$count_data,
                                             pulseData$conditions,
                                             pulseData$spikeins)
  pulseData
}
