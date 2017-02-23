
#' Fit parameters given the initial values and the parameter names
#'
#' @param pd the \code{\link{PulseData}} object
#' @param par the parameter named list
#' @param namesToOptimise a vector of names
#' @param opts a list with optimisation options
#'
#' @return a list with fitted parameters
#' @export
#'
#' @examples
fitParams <- function(pd, par, namesToOptimise, opts) {
  lb <- unlist(opts$lb[namesToOptimise])
  ub <- unlist(opts$ub[namesToOptimise])
  objective <- ll(par = par, namesToOptimise = namesToOptimise, pd = pd)
  x <- unlist(par[namesToOptimise])
  x <- optim(
    x,
    objective,
    method = "L-BFGS-B",
    control = list(parscale = x),
    lower = lb,
    upper = ub,
    counts = pd$counts
  )$par
  relist(x, par[namesToOptimise])
}

fitParamsSeparately <- function(pd, par, namesToOptimise, opts) {
  lb <- unlist(opts$lb[namesToOptimise])
  ub <- unlist(opts$ub[namesToOptimise])
  p <- data.frame(par[namesToOptimise])
  objective <- ll(par, namesToOptimise, pd,  singleValue=TRUE)
  for (i in seq_len(p[, 1])) {
    p[i, ] <- optim(
      p[i, ],
      objective,
      method = "L-BFGS-B",
      control = list(parscale = p[i,]),
      lower = lb,
      upper = ub,
      counts = pd$counts[i,]
    )$par
  }
  as.list(p)
}


#' Fit fraction normalisation coefficients
#'
#' @importFrom  stats optimise
fitNormFactors <- function(pd, par, opts) {
  lb <- unlist(opts$lb$normFactors)
  ub <- unlist(opts$ub$normFactors)
  objective <- llnormFactors(par = par, pd = pd)
  x <- unlist(par$normFactors)
  x <- optim(
    x,
    objective,
    method = "L-BFGS-B",
    control = list(parscale = x),
    lower = lb,
    upper = ub,
    counts = pd$counts
  )$par
  relist(x, par$normFactors)
}


getMaxRelDifference <- function(x,y) max(abs(1 - unlist(x)/unlist(y)))

#' Fit the model by MLE
#'
#' @param pulseData PulseData object
#' @param par initial guess for parameters as a \code{list}.  
#' @param options \code{list} of options
#'
#' @return a list with the fitted parameters in the same form as
#' the initial guess `par`
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
#'    negative binomial distribution, see \code{dnbinom}}
#'    }
#' @export
#'
#' @examples 
#' \dontrun{
#' fitResult <- fitModel(pd, par)
#' }
fitModel <- function(pulseData, par, options) {
  param_names <- names(par$params)
  log2screen(options, cat("\n"))
  rel_err <- Inf
  shared_rel_err <- ifelse(is.null(par$shared), 0, Inf)
  fraction_rel_err <- ifelse(is.null(par$fraction_factors), 0, Inf)
  while (rel_err > options$tolerance$params ||
         shared_rel_err > options$tolerance$shared ||
         fraction_rel_err > options$tolerance$fraction) {
    # Fit shared params
    if (!is.null(par$shared_params)) {
      shared_params <- fitParams(
        pd = pd,
        par = par,
        namesToOptimise = shared_names,
        opts = options
      )
      shared_rel_err <- getMaxRelDifference(shared_params, par[shared_names])
      par$shared <- shared_params
    }
    # Fit params for every genes individually
    params <- fitGeneParameters(pulseData, par, options)
    rel_err <- getMaxRelDifference(params, par$params)
    par$params <- params
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
  par
}
