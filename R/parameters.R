
# default params
# - tolerance
# - boundaries
# - 
.defaultParams <-  list(
  tolerance = list(
    params = 1e-3,
    shared = 1e-2,
    fraction = 1e-3
  ),
  verbose = "silent",
  cores = 1,
  lb = list(size = 10),
  ub = list(size = 1e10)
)

addDefault <- function(options) {
  nonSpecified <- setdiff(names(.defaultParams), names(options))
  options[nonSpecified] <- .defaultParams[nonSpecified]
  options
}

setKnownGeneSpecific <- function(known, par){
  par$known <- known
  par
}

plist <- function(params = NULL,
                  shared = NULL,
                  fraction_factors = NULL,
                  size = NULL) {
  list(
    params = params,
    shared = shared,
    fraction_factors = fraction_factors,
    size = size
  )
}


validateOptions <- function(o){
  if (!is.list(o))
    stop("Options must be a list")
  if (is.null(o$cores) || o$cores < 1)
    stop("Please specify correct number of cores")
  checkThresholds(o)
  checkBoundaries(o)
}

# aligns lb and ub only if they are named
# skip "size" parameter
alignBoundaries <- function(options){
  stopifnot(all(sort(names(options$lb)) == sort(names(options$ub))))
  param_names <- names(options$lb)
  param_names <- param_names[param_names != "size"]
  for (p in param_names) {
    lnames <- names(options$lb[[p]])
    unames <- names(options$ub[[p]])
    if (!is.null(lnames) && !is.null(unames)) {
      options$ub[[p]] <- options$ub[[p]][lnames]
    } else
      if (xor(is.null(lnames), is.null(unames)))
        stop(paste("One of boundaries is not named in ", p))
  }
  options
}


checkBoundaries <- function(options) {
  lnames <- names(options$lb)
  unames <- names(options$ub)
  missingBoundaries <- setdiff(union(lnames, unames), intersect(lnames, unames))
  if (length(missingBoundaries) > 0)
    stop(paste("Please specify missing boundaries for \n", missingBoundaries))
  hasEqualLength <- vapply(lnames,
                           function(p) {
                             length(options$lb[[p]]) == length(options$ub[[p]])
                           }, logical(1))
  if (!(all(hasEqualLength)))
    stop(paste("Length of upper and lower boundaries are not equal for ",
               lnames[!hasEqualLength]))
  correctValues <- vapply(lnames,
                          function(p) {
                            all(options$lb[[p]] < options$ub[[p]])
                          }, logical(1))
  if (!all(correctValues))
    stop(paste( "Lower boundaries values are not less than the upper ones in ",
      lnames[!correctValues]))
  alignBoundaries(options)
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

checkPlist <- function(plist) {
  if (is.null(plist))
    stop("parameter list must be not NULL")
  if (!is.list(plist))
    stop("parameter list must be a list")
  if (length(plist) > 0) {
    listNames <- c("params", "shared", "size", "fraction_factor")
    n <- setdiff(names(plist), listNames)
    if (is.null(n) || length(n) > 0)
      stop(paste(
        "parameter list can have only the following named items: \n",
        paste(listNames, collapse = ", ")
      ))
  }
}

#' Set optimization boundaries for the model parameters.
#'
#' @param params a list of gene-specific parameter boundaries
#' @param shared a list of shared parameters boundaries
#' @param fraction_factors the  lower and the upper boundaries 
#'        for the fraction factors
#' @param size the lower and the upper boundaries for the size parameter for
#' \code{\link{dnbinom}}.
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
setBoundaries <- function(params = NULL,
                          shared = NULL,
                          fraction_factors = NULL,
                          size = NULL,
                          options = .defaultParams) {
  if (!is.list(options))
    stop("Options must be a list")
  options <- addDefault(options)
  plist <- as.list(match.call())[-1]
  plist$options <- NULL
  plist <- lapply(plist, eval)
  .getBoundaryValue <- function(x, i) {
    x <- unlist(x)
    if (length(x) != 2)
      stop("Boundaries must have length of 2")
    x[i]
  }
  for (type in names(plist)) {
    if (is.list(plist[[type]])) {
      # add try
      options$lb[[type]] <- vapply(plist[[type]],
                                   .getBoundaryValue,
                                   i = 1,
                                   FUN.VALUE = double(1))
      options$ub[[type]] <- vapply(plist[[type]],
                                   .getBoundaryValue,
                                   i = 2,
                                   FUN.VALUE = double(1))
    } else {
      if (is.vector(plist[[type]]) && length(plist[[type]])) {
        options$lb[[type]] <- plist[[type]][1]
        options$ub[[type]] <- plist[[type]][2]
      }
    }
  }
  options <- alignBoundaries(options)
  checkBoundaries(options)
  options
}

#' Set the stopping criteria in a form of the relative 
#' changes during fitting iterations.
#'
#' @param params a threshold for gene-specific parameter boundaries
#' @param shared a threshold for shared parameters boundaries
#' @param fraction_factors  a threshold for the fraction factors
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
setTolerance <- function(params = NULL,
                         shared = NULL,
                         fraction_factors = NULL,
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
#' @param pulseData a \code{\link{PulseData}} object 
#' @param options an options object
#' @param params a data.frame
#' @param shared a named list
#' @param fraction_factors a vector
#' @param size a double; size parameter for the \code{\link{dnbinom}}
#' 
#' @importFrom stats runif
#'
#' @return  a list to provide to the function \code{\link{fitModel}}.
#' @export
#'
initParameters <- function(pulseData,
                           options,
                           params = NULL,
                           shared = NULL,
                           fraction_factors = NULL,
                           size = NULL) {
  res <- list()
  validateOptions(options)
  params <- initGeneParams(options, pulseData, params)
  fraction_factors <- initFractions(options, pulseData, fraction_factors)
  shared <- initShared(options, shared)
  size <- runif(1, options$lb$size, options$ub$size)
  par <- plist(
    params = params,
    shared = shared,
    fraction_factors = fraction_factors,
    size = size
  )
  validateNames(par, options)
  stopIfNotInRanges(par, options)
  par
}

validateNames <- function(par, options){
  options$lb$params <- validate(par$params, options$lb$params)
  options$ub$params <- validate(par$params, options$ub$params)
  if (!is.null(par$shared)) {
    options$lb$shared <- validate(par$shared, options$lb$shared)
    options$ub$shared <- validate(par$shared, options$ub$shared)
  }
  options
}

# if not set - sample
initFractions <- function(options, pulseData, fraction_factors) {
  if (is.null(pulseData$fraction))
    return(NULL)
  if (is.null(options$lb$fraction_factors) ||
      is.null(options$ub$fraction_factors)) 
    stop(
      paste(
        "Fractions are set in the PulseData object, ",
        "but boundaries for the fraction factors are not defined"
      )
    )
  fractionNum <- length(levels(pulseData$fraction))
  if (!is.null(fraction_factors)) {
    if (is.vector(fraction_factors) &&
        length(fraction_factors) == 1) {
      fraction_factors <- rep(fraction_factors, fractionNum)
    } else {
      fraction_factors <- fraction_factors
    }
  } else {
      lb <- options$lb$fraction_factors
      ub <- options$ub$fraction_factors
      if (length(lb) == 1)
        lb <- rep(lb, fractionNum)
      if (length(ub) == 1)
        ub <- rep(ub, fractionNum)
      names(ub) <- names(lb) <- levels(pulseData$fraction)
      fraction_factors <- sampleParams(lb, ub, "fraction_factors")
  }
  fraction_factors
}

initShared <- function(options, shared){
  if (!is.null(shared)) {
    if (length(shared) != length(options$lb$shared))
      if (length(shared) == 1)
        shared <- rep(shared, length(options$lb$shared))
      else
        stop("Wrong number of shared parameters specified")
  } else {
    if (!is.null(options$lb$shared) && !is.null(options$ub$shared)) {
      shared <- sampleParams(options$lb$shared, options$ub$shared, "shared") 
    }
  }
  shared
}

initGeneParams <- function(options, pulseData, params){
  geneNum <- dim(pulseData$count_data)[1]
  if (!is.null(params)) {
    correctLength <- sapply(params, length) %in% c(1, geneNum)
    if (!all(correctLength))
      stop(paste(
        "Length of gene-specific parameters is not correct for: ",
        paste(names(params)[!correctLength], collapse = ", ")
      ))
    result <- data.frame(params)
    # correct row number
    if (dim(result)[1] == 1)
      result <- result[rep(1, geneNum),]
  } else {
    result <- replicate(
      geneNum,
      sampleParams(options$lb$params, options$ub$params, "params"),
      simplify = FALSE)
    result <- do.call(rbind, result)
  }
  rownames(result) <- rownames(pulseData$count_data)
  as.data.frame(result)
}

# checks  if parameters are named and
# orders   boundaries appropriately
# if a a boundary is a double, returns b
validate <- function(p, b) {
  if (is.vector(b) && length(b) == 1 && is.null(names(b)))
    return(b)
  if (is.null(names(p)))
    stop("parameters are not named")
  if (length(p) != length(b))
    stop(
      paste(
        "Number of parameters and the number of boundaries set are",
        "not equal:\n parameter:",
        paste0(names(p), collapse = " "),
        "\n boundaries: ",
        paste0(names(b), collapse = " ")
      )
    )
  if (is.null(names(b))) {
    message("Boundaries for the parameters are not named")
    message("The order is derived from the parameter values:")
    message(paste(names(p), collapse = "; "))
    names(b) <- names(p)
  }
  b <- b[names(p)]
  if (anyNA(names(b)))
    stop(paste("Boundaries for ", names(p)[is.na(names(b))],
               "are not set"))
  b
}

#' Validate list of parameters according to allowed value ranges.
#'
#' @param args a list of parametets
#' @param options an options object
#'
#' @return NULL
#'
stopIfNotInRanges <- function(args, options) {
  is.inRange <- function(x, lb, ub) {
    all(vapply(names(x),
               function(par_name) {
                 all(x[[par_name]] >= lb[[par_name]]) &&
                 all(x[[par_name]] <= ub[[par_name]])
               }, logical(1)))
  }
  inRange <- vapply(names(args)[!is.null(args)],
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
