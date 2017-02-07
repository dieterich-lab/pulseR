
#' Create an object for pulse-change count data
#'
#' @param count_data a matrix; column names correspond to sample names
#' @param conditions a data.frame;
#'   the first column corresponds to the conditions given in \code{formulas}.
#'   May columns named as parameters,
#'   used in the formula definitions, e.g. "time" 
#'   if formula = ~condition + time.
#'   
#' @param formulas a list, created by \code{\link{MeanFormulas}}
#' @param spikeins a vector of characters or indexes, optional;
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
  class(e) <- "PulseData"
  samples <- sort(colnames(count_data))
  e$user_conditions <- conditions[samples,, drop = FALSE]
  e$count_data <- as.matrix(count_data[, samples])
  t <- addKnownShared(formulas, e$user_conditions[samples, , drop = FALSE])
  e$conditions <- t$conditions
  e$formulas <- t$formulas
  e$user_formulas <- formulas
  if (!is.null(spikeins) && !is.null(fractions))
    stop(paste("Fractions can not be specified if spike-ins are given.\n"))
  # create e$fraction if spike-ins are not provided
  if (!is.null(spikeins)) {
    e$spikeins <- spikeins
    e$fraction <- NULL
  } else {
    # if no fractions formula provided, use the whole daat.frame from conditions
    if (!is.null(fractions)) {
      columns <- e$user_conditions[, all.vars(fractions), drop = FALSE]
      e$fraction <- factor(apply(columns, 1, paste, collapse = "."))
    } else {
      e$fraction <- factor(apply(e$user_conditions, 1, paste, collapse = "."))
    }
  }
  e$formulas <- t$formulas
  e$norm_factors <- NULL
  normalise(e)
  e
}

#' @export
print.PulseData <- function(x,...){
  cat("PulseData object")
}

#' Calculate normalisation factors for columns in a matrix
#'
#' @param count_data a matrix; columns correspond to samples.
#'
#' @return a vector of doubles with normalisation factors ordered as columns in 
#'   the \code{count_data}
#'
#' @importFrom  stats median
findDeseqFactorsSingle <- function(count_data)
{
  loggeomeans <- rowMeans(log(count_data))
  deseqFactors <-  apply(count_data, 2, function(x) {
    finitePositive <- is.finite(loggeomeans) & x > 0
    if (any(finitePositive))
      res <- exp(median((log(x) - loggeomeans)[finitePositive], na.rm = TRUE))
    else {
      print(count_data[1:6,])
      stop("Can't normalise accross a condition. 
              Too many zero expressed genes. ")
    }
    res  
  })
  deseqFactors
}

#' Calculate normalisation factors
#'
#' @param count_data integer matrix, colnames correspond to samples 
#'   (rownames in \code{conditions})
#' @param conditions factors to split samples for normalisation
#' @return vector of double; normalisation factors in the same order as 
#'   columns in the \code{count_data}
findDeseqFactorsForFractions <- function(count_data, conditions) {
    deseqFactors <- lapply(
      split(colnames(count_data), conditions),
      function(samples) {
        findDeseqFactorsSingle(count_data[, samples, drop = FALSE])
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
normalise <- function(pulseData) {
  if (!is.null(pulseData$spikeins)) {
    pulseData$norm_factors <-
      findDeseqFactorsSingle(pulseData$count_data[pulseData$spikeins, ])
    genes <- setdiff(rownames(pulseData$count_data),pulseData$spikeins)
    pulseData$count_data <- pulseData$count_data[genes, ]
  } else {
    pulseData$norm_factors <- 
      findDeseqFactorsForFractions(pulseData$count_data, pulseData$fraction)
  }
  pulseData
}

# evaluate formulas in the environment of known params from the conditions
# Returns list( [evaluated formulas], [conditions as vector])
addKnownShared <- function(formulas, user_conditions){
  if (dim(as.matrix(user_conditions))[2] == 1)
    return(list(formulas=formulas, conditions=user_conditions[,1]))
  knownParams <-
    which(colnames(user_conditions) %in% unlist(lapply(formulas, all.vars)))
  conditions <- user_conditions[, c(1,knownParams), drop=FALSE]
  interactions <- interaction(conditions,drop = TRUE, lex.order=TRUE)
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
#'     and shared_params (optional). If \code{fractions} is defined,
#'     \code{par$fraction_factors} must be not \code{NULL}
#' @param conditions a condition data.frame
#' @param fractions a factor for splitting the samples into subsets with 
#'   different normalisation coefficients
#' @return matrix of counts with the order of columns as in conditions 
#' @importFrom stats rnbinom
#' @export
#'
generateTestDataFrom <- function(formulas,
                                 par,
                                 conditions,
                                 fractions = NULL){
  if (!is.null(fractions) && is.null(par$fraction_factors)) {
    stop(
      paste(
        "Fraction factors must be specified in par$fraction_factors according\n",
        "to factor levels in the fractions argument"
      )
    )
  }
  t <- addKnownShared(formulas, conditions)
  formulas <- t$formulas
  conditions_known <- data.frame(condition=t$conditions)
  counts <- list()
  for (i in seq_along(par$individual_params[,1])){
    means <- sapply(formulas, eval, 
      c(as.list(par$individual_params[i,]),
        as.list(par$shared_params),
        as.list(par$known[i,, drop=FALSE])))
    # normalise
    norm_factors <- 1
    if (!is.null(fractions)) {
      fraction_indexes <- as.integer(fractions)
      norm_factors <- c(1, par$fraction_factors)[fraction_indexes]
    }
    indexes <- match(conditions_known$condition, names(formulas))
    counts[[i]] <- rnbinom(
      n    = length(conditions$condition),
      mu   = means[indexes] * norm_factors,
      size = par$size)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- rownames(par$individual_params)
  colnames(counts) <- rownames(conditions)
  counts
}
