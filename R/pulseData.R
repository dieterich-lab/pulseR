
#' Create an object for pulse-change count data
#'
#' @param count_data a matrix; column names correspond to sample names
#' @param conditions a data.frame;
#'   the first column corresponds to the conditions given in \code{formulas}.
#'   May also include column "fraction" and columns named as parameters,
#'   used in the formula definitions, e.g. "time".
#'   
#' @param formulas a list, created by \code{\link{MeanFormulas}}
#' @param spikeins a vector of charecters or indexes, optional;
#'  defines which genes to use as a reference for normalisation 
#' @param fractions a formula, e.g. ~ condition + time (if spike-ins are
#' not provided).
#' the names used in the \code{fractions} defines different fractions,
#' which should have distinct coefficients for mean expression fitting.
#'
#' @return an object of class "PulseData"
#' @export
#'
PulseData <- function(count_data,
                      conditions,
                      formulas,
                      spikeins = NULL,
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
  e$norm_factors <- NULL
  class(e) <- "PulseData"
  e
}

#' @export
normalise <- function(x) UseMethod("normalise")

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
  if (is.null(spikeins)) {
    spikeins <- rownames(count_data)
  }
  deseqFactors <- lapply(
    split(colnames(count_data), conditions),
    function(samples) {
      findDeseqFactorsSingle(count_data[spikeins, samples, drop = FALSE])
    })
  names(deseqFactors) <- NULL
  unlist(deseqFactors)[colnames(count_data)]
}

# Performs DESeq normalisation according to first column of *conditions*
# specified by the user (default), or according to *fraction formula*,
# e.g. ~ condition + time
#' Estimate normalisation factors for fraction with the same method as
#' in the DESeq2 package
#'
#' @param pulseData a PulseData object
#'
#' @return the same PulseData object with estimated normalisation factors
#' @export
#'
normalise.PulseData <- function(pulseData) {
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


#' Create a test count data
#'
#' @param formulas a list
#' @param par a list with individual_params(must have), size (must have) 
#'     and shared_params (optional).
#' @param conditions a condition data.frame
#'
#' @return matrix of counts with the order of columns as in conditions 
#' @export
#'
generateTestDataFrom <- function(formulas, par, conditions) {
  counts <- list()
  for(i in seq_along(par$individual_params[,1])){
    means <- sapply(formulas, eval, 
      c(as.list(par$individual_params[i,]),
        as.list(par$shared_params)))
    # normalise
    if(!is.null(conditions$fraction)){
      fraction_indexes <-(as.integer(conditions$fraction))
      means <- means * c(1, par$fraction_factors)[fraction_indexes]
    }
    indexes <- match(conditions$condition, names(formulas))
    counts[[i]] <- rnbinom(
      n    = length(conditions$condition),
      mu   = means[indexes],
      size = par$size)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- rownames(par$individual_params)
  colnames(counts) <- rownames(conditions)
  counts
}
