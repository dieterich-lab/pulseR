
PulseData <- function(count_data,
                      conditions,
                      formulas,
                      spikeins = rownames(count_data),
                      normalisation = NULL) {
  e <- new.env()
  samples <- sort(colnames(count_data))
  e$count_data <- as.matrix(count_data[, samples])
  t <- addKnownShared(formulas, conditions)
  e$conditions <- t$conditions[samples,, drop=FALSE]
  rownames(e$conditions) <- samples
  e$formulas <- t$formulas
  e$spikeins <- spikeins
  e
}

findDeseqFactorsSingle <- function(count_data)
{
  loggeomeans <- rowMeans(log(count_data))
  deseqFactors <-  apply(count_data, 2, function(x) {
    finitePositive <- is.finite(loggeomeans) & x > 0
    if(any(finitePositive))
      res <- exp(median((log(x) - loggeomeans)[finitePositive], na.rm = TRUE))
    else {
      stop("Can't normalise accross a condition. 
              Too many zero expressed genes. ")
    }
    res  
  })
  deseqFactors
}

# conditions  - factor to split samples for normalisation
findDeseqFactors <- function(count_data, conditions, spikeins) {
  deseqFactors <- lapply(split(colnames(count_data), conditions),
                         function(samples) {
                           findDeseqFactorsSingle(count_data[spikeins, samples, drop=FALSE])
                         })
  names(deseqFactors) <- NULL
  unlist(deseqFactors)[colnames(count_data)]
}

normalise <- function(pulseData, fractions) {
  if (missing(fractions)) {
    splitting_factor <- as.data.frame(pulseData$conditions)[,1]
  } else {
    splitting_factor <- pulseData$conditions[, all.vars(fractions), drop=FALSE]
    splitting_factor <- apply(conditions, 1, paste, collapse = ".")
  }
  pulseData$norm_factors <- findDeseqFactors(
    pulseData$count_data,
    splitting_factor,
    pulseData$spikeins)
  pulseData
}

addKnownShared <- function(formulas, conditions){
  if (dim(as.matrix(conditions))[2] == 1)
    return(
      list(formulas = formulas,
        conditions  = conditions))
  knownParams <- which( colnames(conditions) %in% unlist(lapply(formulas, all.vars)))
  conditions <- conditions[, c(1,knownParams)]
  interactions <- interaction(conditions,drop = TRUE)
  names(interactions)   <- rownames(conditions)
  conditions <- unique(conditions)
  evaledFormulas <- lapply(seq_along(conditions[,1]), function(i){
    substitute_q(formulas[[as.character(conditions[i,1])]],
      as.list(conditions[i,-1, drop=FALSE]))
      })
  names(evaledFormulas) <- unique(interactions)
  list(formulas = evaledFormulas,
    conditions  = interactions)
}
