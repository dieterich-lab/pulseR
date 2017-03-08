

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
  opts <- normaliseBoundaries(opts, par, pd)
  # garantee that boundaries are in the same order as the params
  lb <- unlist(opts$lb[namesToOptimise])
  ub <- unlist(opts$ub[namesToOptimise])
  objective <- ll(par = par, namesToOptimise = namesToOptimise, pd = pd)
  x <- unlist(par[namesToOptimise])
  x <- stats::optim(
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
fitParamsSeparately <- function(pd, par, knownNames, namesToOptimise, opts) {
  opts <- normaliseBoundaries(opts, par, pd)
  # garantee that boundaries are in the same order as the params
  lb <- as.data.frame(opts$lb[namesToOptimise])
  ub <- as.data.frame(opts$ub[namesToOptimise])
  p <- data.frame(par[namesToOptimise])
  objective <- ll(par, namesToOptimise, pd, byOne = TRUE)
  par[namesToOptimise] <- NULL
  fixedPars <- par
  for (i in seq_along(p[, 1])) {
    fixedPars[knownNames] <- lapply(par[knownNames], `[[`, i)
    p[i, ] <- stats::optim(
      unlist(p[i, ]),
      objective,
      method = "L-BFGS-B",
      control = list(parscale = p[i,]),
      lower = lb[i,],
      upper = ub[i,],
      counts = pd$counts[i,],
      fixedPars = fixedPars
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
  x <- stats::optim(
    x,
    objective,
    method = "L-BFGS-B",
    control = list(parscale = (lb+ub)/2),
    lower = lb,
    upper = ub,
    counts = pd$counts
  )$par
  relist(c(1,x), par$normFactors)
}

getMaxRelDifference <- function(x, y, eps)
{
  max(abs(1 - unlist(x) / (unlist(y))), na.rm = TRUE)
}

#' Fit the model by MLE
#'
#' @param pulseData PulseData object
#' @param par initial guess for parameters as a \code{list}.  
#' @param options \code{list} of options
#'
#' @return a list with the fitted parameters in the same form as
#' the initial guess `par`
#'     
#' @export
#'
#' @examples 
#' \dontrun{
#' fitResult <- fitModel(pd, par)
#' }
fitModel <- function(pulseData, par, options){
  known <- setdiff(names(par), names(options$lb))
  toFit <- setdiff(names(par), c("size", "normFactors", known))
  options <- normaliseBoundaries(
    options, par[setdiff(names(par), known)], pulseData)
  len <- vapply(par, length, integer(1))
  sharedParams <- toFit[len[toFit] == 1] 
  geneParams <- toFit[len[toFit] > 1]
  knownGenePars <- names(len[known] > 1)
  if (!is.null(pulseData$interSampleCoeffs) && is.null(par$normFactors)) {
    par$normFactors <- assignList(pulseData$interSampleCoeffs, 1)
  }
  log2screen(options, cat("\n"))
  rel_err <- Inf
  shared_rel_err <- ifelse(length(sharedParams) == 0, 0, Inf)
  fraction_rel_err <- ifelse(is.null(par$normFactors), 0, Inf)
  while (rel_err > options$tolerance$params ||
         shared_rel_err > options$tolerance$shared ||
         fraction_rel_err > options$tolerance$normFactors) {
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
      knownNames = knownGenePars,
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
    knownNames = knownGenePars,
    opts = options
  )
  par[geneParams] <- res
  par
}

