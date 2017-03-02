
# default params
# - tolerance
# - boundaries
# - 
.defaultParams <-  list(
  tolerance = list(
    params = 1e-3,
    shared = 1e-2,
    normFactors = 1e-3
  ),
  verbose = "silent",
  cores = 1,
  lb = list(size = 10),
  ub = list(size = 1e10)
)

# extend boundaries to param length
.b <- function(b, par) {
  for (p in names(b)) {
    if (length(b[[p]]) == 1) 
      b[[p]] <- rep(b[[p]], length(unlist(par[[p]])))
      if (is.list(par[[p]])) 
        b[[p]] <- utils::relist(b[[p]], par[[p]])
  }
  b
}

addNormFactorsBoundaries <- function(options, pd){
  options$lb$normFactors <- pulseR:::assignList(pd$interSampleCoeffs, .01)
  options$ub$normFactors <- pulseR:::assignList(pd$interSampleCoeffs, 20)
  options
}

normaliseBoundaries <- function(options, par, pd){
  if (!is.null(pd$interSampleCoeffs))
    options <- addNormFactorsBoundaries(options, pd)
  toExtend <- setdiff(names(options$lb), "normFactors")
  options$lb[toExtend] <- .b(options$lb[toExtend], par[toExtend])
  options$ub[toExtend] <- .b(options$ub[toExtend], par[toExtend])
  options
}

correctLength <- function(options, par){
 options$lb <- .b(options$lb, par)
 options$ub <- .b(options$ub, par)
 options
}

addDefault <- function(options) {
  nonSpecified <- setdiff(names(.defaultParams), names(options))
  options[nonSpecified] <- .defaultParams[nonSpecified]
  options
}

validateOptions <- function(o){
  if (!is.list(o))
    stop("Options must be a list")
  if (is.null(o$cores) || o$cores < 1)
    stop("Please specify correct number of cores")
  checkThresholds(o)
}

checkThresholds <- function(options){
  if (is.null(options$tolerance))
    stop("No tolerance is specified")
  isValid <- vapply(names(options$tolerance),
                    function(p) {
                      if (is.vector(options$tolerance[[p]])   &&
                          length(options$tolerance[[p]]) == 1 &&
                          is.numeric(options$tolerance[[p]]))
                        if (options$tolerance[[p]] > 0)
                          return(TRUE)
                      FALSE
                    }, logical(1))
  if (!all(isValid))
    stop("Tolerance must be a single positive number")
}


#' Set optimization boundaries for the model parameters.
#'
#' @param b
#' @param options an options object to use as a basis for a new parameter set
#'
#' @return   an options object with the new parameter values
#' @details If no options object is provided, the default value is used.
#' Boundaries are provided as a named list of vectors or lists with 
#' the length 2, see example.
#' @export
#' 
#' @examples
#' setBoundaries(params = list(a = c(1,2), b = c(10, 20)))
#'
setBoundaries <- function(b, options = .defaultParams) {
  if (!is.list(options))
    stop("Options must be a list")
  options <- addDefault(options)
  .getBoundaryValue <- function(x, i) {
    x <- unlist(x)
    if (length(x) != 2)
      stop("Boundaries must have length of 2")
    x[i]
  }
  for (p in names(b)) {
        options$lb[[p]] <- b[[p]][1]
        options$ub[[p]] <- b[[p]][2]
  }
  options
}

#' Set the stopping criteria in a form of the relative 
#' changes during fitting iterations.
#'
#' @param params a threshold for gene-specific parameter boundaries
#' @param shared a threshold for shared parameters boundaries
#' @param normFactors a threshold for the fraction factors
#' 
#' @param options an options object to use as a basis for a new parameter set
#'
#' @return   an options object with the new parameter values
#' @details If no options object is provided, the default options are used
#' as a base.
#' A threshold  represents the relative changes in parameter values 
#' between two subsequent fitting iterations.
#' @export
#' @examples 
#' setTolerance(params = 1e-2)
#'
setTolerance <- function(params = .01,
                         shared = .01,
                         normFactors = .01,
                         options = .defaultParams) {
  if (!is.list(options))
    stop("Options must be a list")
  options <- addDefault(options)
  plist <- as.list(match.call())[-1]
  plist$options <- NULL
  plist <- lapply(plist, eval)
  options$tolerance[names(plist)] <- plist
  checkThresholds(options)
  options
}

#' Specify fitting options
#'
#' @param verbose if "verbose" relative changes of parameters 
#' between two fitting iterations are printed
#' @param cores \code{integer}, natural number
#' @param options an option object
#'
#' @return an option object with modified specified parameters.
#' @export
#'
setFittingOptions <- function(
    verbose = c("silent", "verbose"),
    cores = 1,
    options
    ){
  if (missing(options))
    options <- .defaultParams
  if (!missing(verbose))
    match.arg(verbose)
  args <- as.list(match.call())[-1]
  args$options <- NULL
  args <- lapply(args, eval)
  options[names(args)] <- args
  options <- addDefault(options)
  validateOptions(options)
  options
}


sampleParams <- function(lb, ub, paramName) {
  n <- max(length(lb), length(ub))
  result <- runif(n, lb, ub)
  if (is.null(names(lb)) && is.null(names(ub)))
    stop(paste("Please provide names for the boundaries in ", paramName))
  if (!is.null(names(lb)))
    names(result) <- names(lb)
  else
    names(result) <-  names(ub)
  result
}

#' Initialize first guess for the parameters 
#'
#' @param par a list with parameter values
#' @param pulseData a \code{\link{PulseData}} object 
#' @param options an options object
#' 
#' @importFrom stats runif
#'
#' @return  a list to provide to the function \code{\link{fitModel}}.
#' @export
#'
initParameters <- function(par, geneParams, pulseData, options) {
  validateOptions(options)
  nGenes <- dim(pulseData$counts)[1]
  for (g in geneParams) {
    if (is.null(par[[g]])) {
      par[[g]] <-  runif(nGenes, options$lb[[g]], options$ub[[g]])
    } else {
      if (length(par[[g]]) == 1)
        par[[g]] <- rep(par[[g]], nGenes)
    }
  }
  notSet <- setdiff(names(options$lb), c("normFactors", names(par)))
  notSet <- notSet[!notSet %in% geneParams]
  for (p in notSet) {
    par[[p]] <- runif(1, options$lb[[p]], options$ub[[p]])
  }
  if (!is.null(pulseData$interSampleCoeffs)) {
    if (is.null(par$normFactors))
      par$normFactors <- pulseData$interSampleCoeffs
  }
  stopIfNotInRanges(par, options)
  par
}


#' Validate list of parameters according to allowed value ranges.
#'
#' @param args a list of parametets
#' @param options an options object
#'
#' @return NULL
#'
stopIfNotInRanges <- function(args, options) {
  args <- args[names(options$lb)]
  if (!is.null(args$size)) {
    if (args$size < options$lb$size || args$size > options$ub$size)
      stop("Error: Argument 'size' is not within the specified range\n")
    args$size <- runif(1, options$lb$size, options$ub$size)
  }
  options$lb <- .b(options$lb, args)
  options$ub <- .b(options$ub, args)
  is.inRange <- function(x, lb, ub) {
                 all(unlist(x) >= unlist(lb)) &&
                 all(unlist(x) <= unlist(ub))
  }
  inRange <- vapply(names(args),
                    function(p) {
                        is.inRange(args[[p]], options$lb[[p]], options$ub[[p]])
                    }, logical(1))
  if (!all(inRange)) {
    msg <-  sapply(
      names(args)[!inRange],
      function(p) {
        paste0("Error: Argument '", p, "' is not within the specified range\n")
      })
    stop(msg)
  }
}

#' Estimate initial guess for the mean expression level 
#'
#' @param pulseData  the \code{\link{PulseData}} object
#' @param totalLabel a character, the name of the factor level in the
#' condition matrix, which correspond to the total fraction ("total" by default).
#' @param fun a function used to estimate the expression level. 
#' Possible variants are mean, median and adjusted geometric mean (i.e.
#' $exp(mean(log(x + .5)))$.
#' 
#'
#' @return a vector of expression level estimations for every gene
#' 
#' @details Use this function to estimate mean read numbers
#' on the basis of the total fraction.
#' @export
#'
guessMeans <- function(pulseData,
                       totalLabel = "total",
                       fun = c("mean", "geomean", "median")) {
  fun <- match.arg(fun)
  fun <- switch(
    fun,
    mean = mean,
    median = median,
    geomean = function(x)
      exp(mean(log(x + .01)))
  )
  totals <- pulseData$user_conditions[,1] == totalLabel
  apply(pulseData$count_data, 1, fun)
}