

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
                                knownNames,
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
    fixedPars[knownNames] <- lapply(par[knownNames], `[[`, i)
    p[i,] <- .fitGene(p, i, objective, lb, ub, fixedPars, pd$counts)$par
  }
  as.list(p)
}

#' fit params for i-th gene
#' p is a data.frame with the being fitted parameters by column
#' objective is a function to optimise 
#' the calling convention if f(x, counts, fixedPars),
#' where x are the parameters to fit, fixed is a character vector of gene-
#' sepecific parameters which are fixed
.fitGene <- function(p, i, objective, lb, ub, fixedPars, counts) {
  stats::optim(
    unlist(p[i,]),
    objective,
    method = "L-BFGS-B",
    control = list(parscale = p[i, ]),
    lower = lb[i, ],
    upper = ub[i, ],
    counts = counts[i, ],
    fixedPars = fixedPars
  )
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
  relist(c(1,x), par$normFactors)
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
  known <- .getKnownNames(par, options)
  options <- normaliseBoundaries(
    options, par[setdiff(names(par), known)], pulseData)
  toFit <- setdiff(names(par), c("size", "normFactors", known))
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
        options = options
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
      options = options
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
      options = options
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
    if (!is.null(options$resultRDS)) {
      saveRDS(object = par, file = options$resultRDS)
    }
  }
  ## fit gene specific final parameters
  res <- fitParamsSeparately(
    pd = pulseData,
    par = par,
    namesToOptimise = geneParams,
    knownNames = knownGenePars,
    options = options
  )
  par[geneParams] <- res
  if (!is.null(options$resultRDS)) {
    saveRDS(object = par, file = options$resultRDS)
  }
  par
}

# if a paramaeter is not mentioned in the boundaries, 
# it is assumed to be fixed
.getKnownNames <- function(par, options) {
  setdiff(names(par), names(options$lb))
}