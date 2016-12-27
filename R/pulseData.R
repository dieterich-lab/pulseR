
# Create a data structure
# - user_conditions - a df, 1st column corresponds to  initial formulas
# - formulas evaluated using known parameters
# - user_formulas - initial formulas
# - spikein list (default = all genes)
# - conditions - for internal usage - a vector corresponding to  evaluated formulas
# - fraction (default = NULL) - a vector for mapping to fraction_factors
# - count_data - a matrix of counts. colnames<->samples, rownames<->genes
PulseData <- function(count_data,
                      conditions,
                      formulas,
                      spikeins = rownames(count_data),
                      fractions = NULL) {
  e <- new.env()
  samples <- sort(colnames(count_data))
  e$user_conditions <- conditions[samples,,drop=FALSE]
  e$count_data <- as.matrix(count_data[, samples])
  t <- addKnownShared(formulas, e$user_conditions[samples,,drop=FALSE])
  e$conditions <- t$conditions
  e$formulas <- t$formulas
  e$user_formulas <- formulas
  if(!is.null(fractions)){
    columns <- e$user_conditions[, all.vars(fractions), drop=FALSE]
    e$fraction <- apply(columns, 1, paste, collapse = ".")
  }
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

# Performs DESeq normalisation according to first column of *conditions*
# specified by the user (default), or according to *fraction formula*,
# e.g. ~ condition + time
normalise <- function(pulseData) {
  if (is.null(pulseData$fraction)) {
    splitting_factor <- as.data.frame(pulseData$user_conditions)[,1]
  } else {
    splitting_factor <- pulseData$fraction
  }
  pulseData$norm_factors <- findDeseqFactors(
    pulseData$count_data,
    splitting_factor,
    pulseData$spikeins)
  pulseData
}

# evaluate formulas in the environment of known params from the conditions
# Returns list( [evaluated formulas], [conditions as vector])
addKnownShared <- function(formulas, user_conditions){
  if (dim(as.matrix(user_conditions))[2] == 1)
    return(
      list(formulas = formulas,
        conditions  = user_conditions[,1]))
  knownParams <- which( colnames(user_conditions) %in% unlist(lapply(formulas, all.vars)))
  conditions <- user_conditions[, c(1,knownParams)]
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
