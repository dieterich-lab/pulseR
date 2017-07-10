

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
#' @rdname fit
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
#' @rdname fit
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
    p[i,] <- .fitGene(p[i,], i, objective, lb, ub, fixedPars, pd$counts)$par
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
#' @rdname fit
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
#' @rdname fit
#'
#' @examples 
#' \dontrun{
#' fitResult <- fitModel(pd, par)
#' }
fitModel <- function(pulseData, par, options){
  log2screen(options, cat("\n"))
  # identify what to fit and what is fixed
  known   <- .getKnownNames(par, options)
  knownGenePars <- .getKnownGeneNames(par, known) 
  fitSets <- list(
    params = .getGeneToFitNames(par, known),
    shared = .getSharedNames(par, known),
    normFactors = "normFactors"
  )
  if (length(fitSets$shared) == 0)
    fitSets$shared <- NULL
  if (is.null(par$normFactors))
    fitSets$normFactors <- NULL
  # prepare functions and boundaries for optimisation
  options <- normaliseBoundaries(
    options, par[setdiff(names(par), known)], pulseData)
  funs <- list(
    params = function(par) 
      fitParamsSeparately(pulseData, par, knownGenePars, fitSets$params, options),
    shared = function(par)
      fitParams(pulseData, par, fitSets$shared, options),
    normFactors = function(par) 
      list(normFactors = fitNormFactors(pulseData, par, options))
  )
  
  err <- c(params = Inf, shared = Inf, normFactors = Inf)[names(fitSets)]
  while (any(err > unlist(options$tolerance[names(fitSets)]))) {
    for (paramSet in names(fitSets)) {
      parNames <- fitSets[[paramSet]]
      res <- funs[[paramSet]](par)
      err[[paramSet]] <- getMaxRelDifference(res, par[parNames])
      par[parNames] <- res
    }
    par["size"] <- fitParams(pulseData, par, "size", options)
    log2screen(options, progressString(err))
    if (!is.null(options$resultRDS)) {
      saveRDS(object = par, file = options$resultRDS)
    }
  }
  ## fit gene specific final parameters
  par[fitSets$params] <- funs[["params"]](par)
  if (!is.null(options$resultRDS)) {
    saveRDS(object = par, file = options$resultRDS)
  }
  par
}

progressString <- function(err) {
  str <- format(unlist(err),
                digits = 2,
                width = 6)
  paste0("Max Rel.err. in ", 
         paste("[",names(str), str,"]", collapse = " "), "\n")
}

# if a parameter is not mentioned in the boundaries, 
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
