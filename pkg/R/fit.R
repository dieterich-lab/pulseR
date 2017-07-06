

#' Fit parameters given the initial values and the parameter names
#'
#' @param pd the \code{\link{PulseData}} object
#' @param par the parameter named list
#' @param namesToOptimise a vector of names of parameters, which values 
#'   need be optimised
#' @param options a list with optimisation options. For details, see
#' \link{setTolerance}, \link{setFittingOptions}.
#'
#' @return a list with fitted parameters
#' @export
#'
fitParams <- function(pd, par, namesToOptimise, options) {
  options <- normaliseBoundaries(options, par, pd)
  # garantee that boundaries are in the same order as the params
  lb <- unlist(options$lb[namesToOptimise])
  ub <- unlist(options$ub[namesToOptimise])
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
#' The same as \link{fitParams}, but performs optimisation for gene-specific
#' parameters only. Every set of parameters is fitted individually for
#' every gene.
#' 
#' @inheritParams fitParams
#' @param knownNames a vectors of names of the gene-specific parameters, which 
#' are assumed to be fixed during optimisation.
#' @param indexes indexes of genes to fit. By default includes all the genes.
#' @return a list with fitted parameters
#' @export
#'
fitParamsSeparately <- function(pd,
                                par,
                                knownGenePars,
                                namesToOptimise,
                                options, 
                                indexes = seq_len(dim(pd$counts)[1])) {
  options <- normaliseBoundaries(options, par, pd)
  # garantee that boundaries are in the same order as the params
  lb <- as.data.frame(options$lb[namesToOptimise])
  ub <- as.data.frame(options$ub[namesToOptimise])
  p <- data.frame(par[namesToOptimise])
  objective <- ll(par, namesToOptimise, pd, byOne = TRUE)
  par[namesToOptimise] <- NULL
  fixedPars <- par
  for (i in indexes) {
    fixedPars[knownGenePars] <- lapply(par[knownGenePars], `[[`, i)
    p[i,] <- .fitGene(p, i, objective, lb, ub, fixedPars, pd$counts)$par
  }
  as.list(p)
}


#' Fit fraction normalisation coefficients
#' 
#' @inheritParams fitParams
#' 
#' @importFrom  stats optimise
#' @return a list of normalisation factors
#' @export
#' 
fitNormFactors <- function(pd, par, options) {
  lb <- unlist(options$lb$normFactors)[-1]
  ub <- unlist(options$ub$normFactors)[-1]
  objective <- llnormFactors(par = par, pd = pd)
  x <- unlist(par$normFactors)[-1]
  x <- stats::optim(
    x,
    objective,
    method = "L-BFGS-B",
    control = list(parscale = (lb + ub) / 2), 
    lower = lb,
    upper = ub,
    counts = pd$counts
  )$par
  result <- relist(c(1,x), par$normFactors)
  names(result) <- names(pd$interSampleCoeffs)
  result
}

getMaxRelDifference <- function(x, y)
{
  max(abs(1 - unlist(x) / (unlist(y))), na.rm = TRUE)
}

#' Fit the model by MLE
#'
#' @param pulseData a \link{PulseData} object
#' @param par a list with an initial parameters values.
#'   Gene-specific parameters must be set as vectors of the length equal to 
#'   the gene number. \link{initParameters} may simplify the process of initial
#'   values randomisation.
#'   
#'   The list must include `size` parameter for the negative binomial 
#'   distribution and a list with the normalisation factors, if 
#'   no spike-ins are used in the experiment.
#'   
#' @param options a list of options. For more details, see \link{setBoundaries},
#'   \link{setTolerance}, \link{setFittingOptions}
#'
#' @return a list with the fitted parameters in the same form as
#' the initial guess `par` list
#'     
#' @export
#'
#' @examples 
#' \dontrun{
#' fitResult <- fitModel(pd, par)
#' }
fitModel <- function(pulseData, par, options){
  known   <- .getKnownNames(par, options)
  options <- normaliseBoundaries(
    options, par[setdiff(names(par), known)], pulseData)
  sharedParams  <- .getSharedNames(par, known) 
  geneParsToFit <- .getGeneToFitNames(par, known) 
  knownGenePars <- .getKnownGeneNames(par, known) 
  if (!is.null(pulseData$interSampleCoeffs) && is.null(par$normFactors)) {
    par$normFactors <- assignList(pulseData$interSampleCoeffs, 1)
  }
  log2screen(options, cat("\n"))
  funs <- list(
    params = function(par) 
      fitParamsSeparately(pulseData, par, knownGenePars,  geneParsToFit, options),
    shared = function(par)
      fitParams(pulseData, par, sharedParams, options),
    normFactors = function(par) 
      list(normFactors = fitNormFactors(pulseData, par, options))
  )
  sets <- list(
    params = geneParsToFit, shared = sharedParams, normFactors = "normFactors")
  if (length(sharedParams) == 0)
    sets$shared <- NULL
  if (is.null(par$normFactors))
    sets$normFactors <- NULL
  
  err <- c(params = Inf, shared = Inf, normFactors = Inf)
  while (any(err[names(sets)] > unlist(options$tolerance[names(sets)]))) {
    # Fit shared params
    for (paramSet in names(sets)) {
      parNames <- sets[[paramSet]]
      res <- funs[[paramSet]](par)
      err[[paramSet]] <- getMaxRelDifference(res, par[parNames])
      par[parNames] <- res
    }
    par["size"] <- fitParams(pulseData, par, "size", options)
    log2screen(options,progressString(err))
    if (!is.null(options$resultRDS)) {
      saveRDS(object = par, file = options$resultRDS)
    }
  }
  ## fit gene specific final parameters
  res <- funs[["params"]](par)
  par[geneParsToFit] <- res
  if (!is.null(options$resultRDS)) {
    saveRDS(object = par, file = options$resultRDS)
  }
  par
}

progressString <- function(err) {
  str <- format(unlist(err),
                digits = 2,
                width = 6)
  paste0("Max Rel.err. in [params: ",
         str[1],
         "]  [shared: ",
         str[2],
         "]  [fractions: ",
         str[3],
         "]    \r")
}

# if a paramaeter is not mentioned in the boundaries, 
# it is assumed to be fixed
.getKnownNames <- function(par, options) {
  setdiff(names(par), names(options$lb))
}

# return names of parameters which must be fitted
.namesToFit <- function(par, known) {
  setdiff(names(par), c("size", "normFactors", known))
}

# return names of gene-specific parameters which must be fitted
# a parameter is assumed to be gene-specific, if its length in `par` is > 1
.getGeneToFitNames <- function(par, known) {
  toFit <- .namesToFit(par, known)
  len <- vapply(par, length, integer(1))
  geneParams <- toFit[len[toFit] > 1]
  geneParams
}

# return names of known gene-specific parameters
# a parameter is assumed to be gene-specific, if its length in `par` is > 1
.getKnownGeneNames <- function(par, known) {
  len <- vapply(par, length, integer(1))
  knownGenePars <- names(len[known] > 1)
  knownGenePars
}

# return names of shared parameters
# a parameter is assumed to be shared, if its length in `par` is 1
.getSharedNames <- function(par, known) {
  toFit <- .namesToFit(par, known)
  len <- vapply(par, length, integer(1))
  sharedParams <- toFit[len[toFit] == 1] 
  sharedParams
}
