
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

validateOptions <- function(o){
  if (is.null(o$cores) || o$cores < 1)
    stop("Please specify correct number of cores")
  if (!all(vapply(o$tolerance, function(x) x > 0, logical(1))))
    stop("Tolerance must be a positive number")
  missingBoundaries <- setdiff(names(o$lb), names(o$ub))
  if (length(missingBoundaries) > 0)
    stop(paste("Please specify missing boundaries for \n", missingBoundaries))
  hasEqualLength <- vapply(names(o$lb), 
         function(p){
           length(o$lb[[p]]) == length(o$ub[[p]])
         }, logical(1))
  if (!(all(hasEqualLength)))
    stop(paste("Length of upper and lower boundaries are not equal for ",
         names(o$lb)[!hasEqualLength]))
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
                          options = NULL) {
  if (missing(options))
    options <- .defaultParams
  validateOptions(options)
  args <- as.list(match.call())[-1]
  args <- lapply(args, eval)
  options$lb[names(args)] <- lapply(args, `[[`, 1)
  options$ub[names(args)] <- lapply(args, `[[`, 2)
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
                         options = NULL) {
  if (missing(options))
    options <- .defaultParams
  args <- as.list(match.call())[-1]
  args <- lapply(args, eval)
  options$tolerance[names(args)] <- args
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
  args <- as.list(match.call())[-1]
  args <- lapply(args, eval)
  options[names(args)] <- args
  validateOptions(options)
  options
}

#' Initialize first guess for the parameters 
#'
#' @param options an options object
#' @param params a data.frame
#' @param shared a named list
#' @param fraction_factors a vector
#'
#' @return
#' a list to provide to the function \code{\link{fitModel}}.
#' @export
#'
initParams <- function(options,
                       params,
                       shared = NULL,
                       fraction_factors = NULL) {
  validateOptions(options)
  args <- as.list(match.call())[-1]
  args$options <- NULL
  args <- lapply(args, eval)
  stopIfNotInRanges(args, options)
  notSpecified <- setdiff(names(args), names(options$lb))
  args
}


#' Validate list of parameters according to allowed value ranges.
#'
#' @param args a list of parametets
#' @param options an options object
#'
#' @return NULL
#'
stopIfNotInRanges <- function(args, options) {
  inRange <- vapply(names(args),
                    function(p) {
                      (args[[p]] > options$lb[[p]]) &&
                      (args[[p]] < options$ub[[p]])
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
 

#' Fit gene-specific parameters
#'
#' @param pulseData PulseData object
#' @param par all needed parameters as a \code{list}
#' @param options \code{list} of options
#'
#' @return fitted parameters as a data.frame ordered as initial parameters
#' @importFrom  stats optim dnbinom
#' @importFrom parallel mclapply
#' 
fitGeneParameters <- function(pulseData, par, options) {
  param_names <- colnames(par$params)
  objective <- ll_gene(pulseData, par)
  new_params <- list()
  new_params <- mclapply(
    X = seq_len(dim(par$params)[1]),
    FUN = function(i) {
      olds <- par$params[i,,drop=FALSE]
      optim(
        olds,
        objective,
        method = "L-BFGS-B",
        lower = options$lb$params,
        upper = options$ub$params,
        control = list(parscale = olds),
        counts = pulseData$count_data[i,],
        known=par$known[i,, drop=FALSE]
      )$par
    }
    ,mc.cores = options$cores
  )
  new_params <- do.call(rbind, new_params)
  rownames(new_params) <- rownames(par$params)
  as.data.frame(new_params)
}

fitSharedParameters <- function(pulseData, par, options) {
  shared_objective <- ll_shared_params(pulseData, par)
  shared_params <- optim(
    unlist(par$shared),
    shared_objective,
    method = "L-BFGS-B",
    lower = options$lb$shared,
    upper = options$ub$shared
  )$par
  names(shared_params) <- names(par$shared)
  as.list(shared_params)
}

#' Fit the dispersion parameter
#'
#' @param pulseData PulseData object
#' @param par all needed parameters as a \code{list}
#' @param options \code{list} of options
#'
#' @return the MLE of the \code{size} for \code{\link{dnbinom}}
#' @importFrom  stats optimise
#'
fitDispersion <- function(pulseData, par, options) {
  dispersion_objective <- ll_dispersion(pulseData, par)
  interval <- c(options$lb$size, options$ub$size)
  size <- optimise(dispersion_objective,
                   interval = interval)$minimum
  size
}

#' Fit fraction normalisation coefficients
#'
#' @inheritParams fitDispersion
#' @return vector of normalisation factors; \code{c(1,norm_factors)}
#'   corresponds to the fractions in \code{pulseData$fractions}
#'
#' @importFrom  stats optimise
fitFractions <- function(pulseData, par, options){
  objective <- ll_norm_factors(pulseData, par)
  fraction_factors <- optim(
    unlist(par$fraction_factors)[-1],
    objective,
    method = "L-BFGS-B",
    lower  = options$lb$fraction_factors,
    upper  = options$ub$fraction_factors
  )$par
  c(1,fraction_factors)
}

getMaxRelDifference <- function(x,y) max(abs(1 - unlist(x)/unlist(y)))

#' Fit the model by MLE
#'
#' @param pulseData PulseData object
#' @param par initial guess for parameters as a \code{list}.  
#' @param options \code{list} of options
#'
#' @return a list `l` with fitted parameters `l$par` and
#'     formulas `l$formulas` 
#'     
#' @details 
#'    In the initial guess \code{par},
#'    the following list items need to be provided:
#'    \itemize{
#'    \item{params: }{list of gene-specific parameters}
#'    \item{shared: }{list of shared parameters}
#'    \item{fraction_factors (if relevant): }{a vector of fraction factors
#'    in the order of \code{levels(pulseDatad$fractions)}} 
#'    \item{size (optional): }{an initial guess for the size parameter of the
#'    negative binomial distribution, see \code{\link{dnbinom}}}
#'    }
#' @export
#'
#' @examples 
#' \dontrun{
#' fitResult <- fitModel(pd, par)
#' }
fitModel <- function(pulseData, par, options = .defaultParams) {
  param_names <- names(par$params)
  log2screen(options, cat("\n"))
  rel_err <- Inf
  shared_rel_err <- ifelse(is.null(par$shared), 0, Inf)
  fraction_rel_err <- ifelse(is.null(par$fraction_factors), 0, Inf)
  while (rel_err > options$tolerance$params||
         shared_rel_err > options$tolerance$shared||
         fraction_rel_err > options$tolerance$fraction) {
    # Fit shared params
    if (!is.null(par$shared_params)) {
      shared_params <- fitSharedParameters(pulseData, par, options)
      shared_rel_err <- getMaxRelDifference(shared_params, par$shared)
      par$shared <- shared_params
    }
    # Fit params for every genes individually
    params <- fitGeneParameters(pulseData, par, options)
    rel_err <- getMaxRelDifference(params, par$params)
    par$params<- params
    if (!is.null(par$fraction_factors)) {
      res <- fitFractions(pulseData, par, options)
      fraction_rel_err <- getMaxRelDifference(res, par$fraction_factors)
      par$fraction_factors <- res
    }
    par$size <- fitDispersion(pulseData, par, options)
    str <- format(c(rel_err, shared_rel_err, fraction_rel_err),
                  digits = 2,
                  width = 6)
    log2screen(options, cat(
      paste0(
        "Max Rel.err. in [params: ", str[1],
        "]  [shared: ", str[2], 
        "]  [fractions: ", str[3],
        "]    \r"
      )
    ))
  }
  ## fit gene specific final parameters
  par$individual_params <- fitGeneParameters(pulseData, par, options)
  names(par$fraction_factors)  <- levels(pulseData$fraction)
  list(par = par, formulas = pulseData$formulas)
}
