
PulseData <- function(count_data,
                      conditions,
                      formulas,
                      spikeins = rownames(count_data)) {
  e <- new.env()
  samples <- sort(colnames(count_data))
  e$count_data <- as.matrix(count_data[, samples])
  e$conditions <- conditions[samples,,drop=FALSE]
  e$formulas <- formulas
  e$spikeins <- spikeins
  e
}

findDeseqFactorsSingle <- function(count_data)
{
  loggeomeans <- rowMeans(log(count_data))
  deseqFactors <-  apply(count_data, 2, function(x) {
    exp(median((log(x) - loggeomeans)[is.finite(loggeomeans) & x > 0],
               na.rm = TRUE))
  })
  deseqFactors
}

# conditions  - factor to split samples for normalisation
findDeseqFactors <- function(count_data,
                             conditions,
                             spikeins) {
  deseqFactors <- lapply(split(names(conditions), conditions),
                         function(samples) {
                           findDeseqFactorsSingle(count_data[spikeins, samples, drop=FALSE])
                         })
  names(deseqFactors) <- NULL
  unlist(deseqFactors)[colnames(count_data)]
}

normalise <- function(pulseData, fractions) {
  if (missing(fractions)) {
    conditions <- pulseData$conditions[, 1]
  } else {
    conditions <- pulseData$conditions[, all.vars(fractions), drop=FALSE]
    conditions <- apply(conditions, 1, paste, collapse = ".")
  }
  names(conditions) <- rownames(pulseData$conditions)
  pulseData$norm_factors <- findDeseqFactors(pulseData$count_data,
                                             conditions,
                                             pulseData$spikeins)
  pulseData
}

addKnownShared <- function(formulas, conditions){
  if (dim(conditions)[2] == 1)
    return(list(formulas = formulas,
                conditions = conditions))
  interactions <- interaction(conditions,drop = TRUE)
  conditions <- unique(conditions)
  evaledFormulas <- lapply(seq_along(conditions[,1]), function(i){
    substitute_q(formulas[[conditions[i,1]]],
                 as.list(conditions[i,-1, drop=FALSE]))
  })
  names(evaledFormulas) <- interactions
  list(formulas = evaledFormulas,
       conditions = interactions)
}
