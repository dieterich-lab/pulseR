
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

#' Fit parameters with separate likelihood functions
#'
#' @inheritParams fitParams
#' @return a list with fitted parameters
#' @export
#'
fitParamsSeparately <- function(pd, par, namesToOptimise, opts) {
  lb <- as.data.frame(opts$lb[namesToOptimise])
  ub <- as.data.frame(opts$ub[namesToOptimise])
  p <- data.frame(par[namesToOptimise])
  objective <- ll(par, namesToOptimise, pd,  singleValue=TRUE)
  for (i in seq_along(p[, 1])) {
    p[i, ] <- optim(
      unlist(p[i, ]),
      objective,
      method = "L-BFGS-B",
      control = list(parscale = p[i,]),
      lower = lb[i,],
      upper = ub[i,],
      counts = pd$counts[i,]
    )$par
  }
  as.list(p)
}


#' Fit fraction normalisation coefficients
#'
#' @importFrom  stats optimise
fitNormFactors <- function(pd, par, opts) {
  lb <- unlist(opts$lb$normFactors)[-1]
  ub <- unlist(opts$ub$normFactors)[-1]
  objective <- llnormFactors(par = par, pd = pd)
  x <- unlist(par$normFactors)[-1]
  x <- optim(
    x,
    objective,
    method = "L-BFGS-B",
    control = list(parscale = x),
    lower = lb,
    upper = ub,
    counts = pd$counts
  )$par
  relist(c(1,x), par$normFactors)
}


fitAll <- function(pd, par, opts) {
  par$normFactors <- par$normFactors[-1]
  opts$lb$normFactors <- opts$lb$normFactors[-1]
  opts$ub$normFactors <- opts$ub$normFactors[-1]
  lb <- unlist(opts$lb[names(par)])
  ub <- unlist(opts$ub[names(par)])
  objective <- totalll(par, pd)
  x <- unlist(par)
  x <- optim(
    x,
    objective,
    method = "L-BFGS-B",
    control = list(parscale = x),
    lower = lb,
    upper = ub,
    counts = pd$counts
  )$par
  x <- relist(x, par)
  x$normFactors <- c(1,x$normFactors)
  x
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
  len <- vapply(par, length, integer(1))
  sharedParams <- names(par)[len == 1 & names(par) != "size"]
  geneParams <- names(par)[len > 1 & names(par) != "normFactors"]
  log2screen(options, cat("\n"))
  rel_err <- Inf
  shared_rel_err <- ifelse(length(sharedParams) == 0, 0, Inf)
  fraction_rel_err <- ifelse(is.null(par$normFactors), 0, Inf)
  while (rel_err > options$tolerance$params ||
         shared_rel_err > options$tolerance$shared ||
         fraction_rel_err > options$tolerance$fraction_factors) {
    # Fit shared params
    if (length(sharedParams) > 0) {
      res <- fitParams(
        pd = pulseData,
        par = par,
        namesToOptimise = sharedParams,
        opts = options
      )
      shared_rel_err <- getMaxRelDifference(res, par[sharedParams])
      par[sharedParams] <- res
    }
    # Fit params for every genes individually
    res <- fitParamsSeparately(
      pd = pulseData,
      par = par,
      namesToOptimise = geneParams,
      opts = options
    )
    rel_err <- getMaxRelDifference(res, par[geneParams])
    par[geneParams] <- res
    if (!is.null(par$normFactors)) {
      res <- fitNormFactors(pulseData, par, options)
      fraction_rel_err <- getMaxRelDifference(res, par$normFactors)
      par$normFactors <- res
    }
    par["size"] <- fitParams(
      pd = pulseData,
      par = par,
      namesToOptimise = "size",
      opts = options
    )
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
  res <- fitParamsSeparately(
    pd = pulseData,
    par = par,
    namesToOptimise = geneParams,
    opts = options
  )
  par[geneParams] <- res
  par
}

