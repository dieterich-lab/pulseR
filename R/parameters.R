
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

plist <- function(params = NULL,
                  shared = NULL,
                  fraction_factors = NULL,
                  size = NULL) {
  as.list(match.call())[-1]
}


validateOptions <- function(o){
  if (!is.list(o))
    stop("Options must be a list")
  if (is.null(o$cores) || o$cores < 1)
    stop("Please specify correct number of cores")
  checkThresholds(o)
  checkBoundaries(o)
  # order of boundaries
}

# aligns lb and ub only if they are named
# skip "size" parameter
alignBoundaries <- function(options){
  stopifnot(sort(names(options$lb) == sort(names(options$ub))))
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


checkBoundaries <- function(options){
  missingBoundaries <- setdiff(names(options$lb), names(options$ub))
  if (length(missingBoundaries) > 0)
    stop(paste("Please specify missing boundaries for \n", missingBoundaries))
  hasEqualLength <- vapply(names(options$lb), 
         function(p){
           length(options$lb[[p]]) == length(options$ub[[p]])
         }, logical(1))
  if (!(all(hasEqualLength)))
    stop(paste("Length of upper and lower boundaries are not equal for ",
         names(options$lb)[!hasEqualLength]))
  
}

checkThresholds <- function(options){
  if (is.null(options$tolerance))
    stop("No tolerance is specified")
  isValid <- vapply(names(options$tolerance),
         function(p) {
           is.vector(options$tolerance[[p]])   &&
           length(options$tolerance[[p]]) == 1 &&
           is.numeric(options$tolerance[[p]])
         }, logical(1))
  if (!all(isValid))
    stop("Tolerance must be a single number")
}

#' Set optimization boundaries for the model parameters.
#'
#' @param params boundaries for gene-specific parameters
#' @param shared boundaries for shared parameters
#' @param fraction_factors boundaries for fraction factors if relevant
#' @param size boundaries for the size parameter of the negative binomial
#' ditribution 
#' @param options an options object to use as a basis for a new parameter set
#'
#' @return   an options object with the new parameter values
#' @details If no options object is provided, the default value is used.
#' @export
#'
setBoundaries <- function(params = NULL,
                          shared = NULL,
                          fraction_factors = NULL,
                          size = NULL,
                          options = .defaultParams) {
  if (!is.list(options))
    stop("Options must be a list")
  options <- addDefault(options)
  args <- as.list(match.call())[-1]
  args <- lapply(args, eval)
  args$options <- NULL
  options$lb[names(args)] <- lapply(args, `[[`, 1)
  options$ub[names(args)] <- lapply(args, `[[`, 2)
  validateOptions(options)
  options
}

#' Set the stopping criteria in a form of the relative 
#' changes during fitting iterations.
#'
#' @param params double
#' @param shared double
#' @param fraction_factors double
#' @param size double
#' @param options an options object to use as a basis for a new parameter set
#'
#' @return   an options object with the new parameter values
#' @details If no options object is provided, the default value is used.
#' @export
#'
setTolerance <- function(params = NULL,
                         shared = NULL,
                         fraction_factors = NULL,
                         size = NULL,
                         options = .defaultParams) {
  if (!is.list(options))
    stop("Options must be a list")
  args <- as.list(match.call())[-1]
  args <- lapply(args, eval)
  args$options <- NULL
  options$tolerance[names(args)] <- args
  isValid <- vapply(names(options$tolerance),
         function(p) {
           is.vector(options$tolerance[[p]])   &&
           length(options$tolerance[[p]]) == 1 &&
           is.numeric(options$tolerance[[p]])
         }, logical(1))
  if (!all(isValid))
    stop("Tolerance must be a single number")
  options <- addDefault(options)
  validateOptions(options)
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
fittingOptions <- function(
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


sampleParams <- function(options, params_name) {
  n <- length(options$lb[[params_name]])
  if (is.null(names(options$lb[[params_name]])))
    stop("Please specify parameter names in boundaries")
  lb <- options$lb[[params_name]]
  ub <- orderBoundaries(names(lb), options$ub[[params_name]])
  result <- runif(n, lb, ub)
  names(result) <- names(lb)
  result
}

#' Initialize first guess for the parameters 
#'
#' @param pulseData a \code{\link{PulseData}} object 
#' @param options an options object
#' @param params a data.frame
#' @param shared a named list
#' @param fraction_factors a vector
#' 
#' @importFrom stats runif
#'
#' @return  a list to provide to the function \code{\link{fitModel}}.
#' @export
#'
initParams <- function(pulseData,
                       options,
                       params = NULL,
                       shared = NULL,
                       fraction_factors = NULL) {
  validateOptions(options)
  args <- as.list(match.call())[-1]
  args$options <- NULL
  args$pulseData <- NULL
  args <- lapply(args, eval, envir = parent.frame()) 
  stopIfNotInRanges(args, options)
  notSpecified <- setdiff(names(options$lb),names(args))
  guess <- lapply(
    notSpecified,
    function(p) {
      if (p != "params") {
        sampleParams(options, p)
      } else {
        geneNum <- dim(pulseData$count_data)[1]
        result <- replicate(geneNum,
                            sampleParams(options, "params"),
                            simplify = FALSE)
        do.call(rbind, result)
      }
    })
  names(guess) <- notSpecified
  args[names(guess)] <- guess
  args
}

# check if names if parameters are named and
# order boundaries appropriately
validate <- function(p, b) {
  if (is.null(names(p)))
    stop("parameters are not named")
  if (is.null(names(b))) {
    message("Boundaries for the parameters are not named")
    message("The order is derived from the parameter values:")
    message(paste(names(p), collapse = "; "))
    names(b) <- names(p)
  }
  b[names(p)]
}

validateNames <- function(args, options){
  for (p in names(args)) {
    options$lb[[p]] <- validate(args[[p]], options$lb[[p]])
    options$ub[[p]] <- validate(args[[p]], options$ub[[p]])
  }
  options
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
                 all(x[[par_name]] > lb[[par_name]]) &&
                 all(x[[par_name]] < ub[[par_name]])
               }, logical(1)))
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
  NULL
}
 
# Align boundary to parameters
# if boundaries are not named - leave them unchanged
.orderBoundaries <- function(plist, b) {
  for (p in names(plist)) {
    if (is.null(b[[p]]))
      stop(paste("Boundaries are not set for the item: ", p))
    if (p != "size" && !is.null(names(b[[p]]))) {
      b <- b[parameter_names]
      if (anyNA(names(b)))
        stop(paste("Boundaries for ", parameter_names[is.na(names(b))],
                   "are not set"))
    }
  }
  b
}

orderBoundaries <- function(plist, options){
 options$lb <- .orderBoundaries(plist, options$lb) 
 options$ub <- .orderBoundaries(puist, options$ub) 
 options
