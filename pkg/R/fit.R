

#' Fit model parameters
#'
#' `fitModel` is the main fitting function which wraps around the others and
#' fit iteratively the whole model. \cr
#' `fitParams` fits only the listed parameters and `fitParamsSeparately` fit
#' gene-specific parameters having the shared parameters (e.g. normalisation
#' factors). \cr
#' `fitNormFactors` fits the normalisation factors having fixed all other
#' parameters.
#'
#' @param pulseData the \code{\link{PulseData}} object
#' @param par a list with an initial parameters values.
#' Names correspond to the parameter names used in formulas.
#' `size` corresponds to the size parameter, `normFactors` stands for the
#' list with normalistion factors (in spike-in free design, see details).
#' There are following parameter types:
#'   - gene-specific parameters must be set as vectors of the length equal to
#'     the gene number. The function \link{initParameters} may
#'     simplify the process of initial values randomisation.
#'   - shared parameters. These are prersented by single numeric values,
#'     which are assumed to be equal between all genes, e.g. additional
#'     normalisation fator.
#'   - the size parameter for the negative binomial
#'     distribution and a list with the normalisation factors, if
#'     no spike-ins are used in the experiment.
#'   - normalisation factors (in case of spike-in free design).
#'
#' @param namesToOptimise a vector of names of parameters, which values
#'   need be optimised
#' @param knownGenePars a vectors of names of the gene-specific parameters,
#'  which  are assumed to be fixed during optimisation.
#' @param indexes indexes of genes to fit. By default includes all the genes.
#' @param options a list of options. For more details, see \link{setBoundaries},
#'   \link{setTolerance}, \link{setFittingOptions}
#'
#' @details If no spike-ins are used, relations between samples are inferred
#' during the model fitting. In this case, the initial parameter list must
#' containg  a field named `normFactors`. The normalistion factors are
#' accepted as a named list, e.g.
#' ```
#' par$normFactors <- list(total_fraction = 1,
#'       pull_down.4 = c(1, 0.01),
#'       pull_down.8 = c(1, 0.01))
#' ```
#' This will define the initial values for the normalisation factors.
#' The very first value is **always** equal 1 irregardless of the user input.
#' This has to be done because the normalisation factors are known
#' only up to some scaling coefficient, because they appear in
#' a multiplication with the expression level or synthesis rate.
#'
#' The structure of the `normFactors` list is identical to the
#' `pulseData$interSampleCoeffs`. This structure is defined by the
#' `formulaIndexes` and `conditions` argumenta in the `PulseData`,
#'  see `\link{PulseData}` for more.
#'
#' `\link{fitParamsSeparately}` is same as \link{fitParams},
#' but performs optimisation for gene-specific parameters only.
#' Every set of parameters is fitted individually for every gene.
#'
#' @return a list with fitted parameters (only which were optimized)
#' @export
#' @rdname fit
#'
fitParams <- function(pulseData, par, namesToOptimise, options) {
  options <- normaliseBoundaries(options, par, pulseData)
  # garantee that boundaries are in the same order as the params
  lb <- unlist(options$lb[namesToOptimise])
  ub <- unlist(options$ub[namesToOptimise])
  parscale <- .5 * (abs(ub) + abs(lb))
  objective <- ll(par = par, namesToOptimise = namesToOptimise, pulseData = pulseData)
  x <- unlist(par[namesToOptimise])
  x <- stats::optim(
    x,
    objective,
    method = "L-BFGS-B",
    control = list(parscale = parscale),
    lower = lb,
    upper = ub,
    counts = pulseData$counts
  )$par
  utils::relist(x, par[namesToOptimise])
}

#' @rdname fit
#' @export
fitParamsSeparately <- function(pulseData,
                                par,
                                knownGenePars,
                                namesToOptimise,
                                options,
                                indexes = seq_len(dim(pulseData$counts)[1])) {
  if (missing(knownGenePars))
    knownGenePars <- character(0)
  if (is.null(options$replicates))
    options$replicates <- 1
  options <- normaliseBoundaries(options, par, pulseData)
  # garantee that boundaries are in the same order as the params
  lb <- as.data.frame(options$lb[namesToOptimise])
  ub <- as.data.frame(options$ub[namesToOptimise])
  p <- data.frame(par[namesToOptimise])
  objective <- ll(par, namesToOptimise, pulseData, byOne = TRUE)
  par[namesToOptimise] <- NULL
  fixedPars <- par
  if (is.null(options$cores))
    options$cores <- 1
  res <- parallel::mclapply(indexes, function(i) {
    fixedPars[knownGenePars] <- lapply(par[knownGenePars], `[[`, i)
    .fitGene(p[i,], i, objective, lb, ub, fixedPars, pulseData$counts,
                      N = options$replicates)$par
  }, mc.cores = options$cores)
  res <- do.call(rbind, res)
  p[indexes,] <- res
  as.list(p)
}


#' @rdname fit
#' @export
fitNormFactors <- function(pulseData, par, options) {
  lb <- unlist(options$lb$normFactors)[-1]
  ub <- unlist(options$ub$normFactors)[-1]
  objective <- llnormFactors(par = par, pulseData = pulseData)
  x <- unlist(par$normFactors)[-1]
  x <- stats::optim(
    x,
    objective,
    method = "L-BFGS-B",
    lower = lb,
    upper = ub,
    counts = pulseData$counts
  )$par
  result <- utils::relist(c(1,x), par$normFactors)
  names(result) <- names(pulseData$interSampleCoeffs)
  result
}

getMaxRelDifference <- function(x, y)
{
  max(abs(1 - unlist(x) / (unlist(y))), na.rm = TRUE)
}


getMaxAbsDifference <- function(x, y)
{
  max(abs(unlist(x) - (unlist(y))), na.rm = TRUE)
}

#' @rdname fit
#' @export
fitModel <- function(pulseData, par, options){
  options <- addDefault(options)
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
  if (is.null(par$normFactors) || options$fixedNorms)
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
  err <- c(err, logLik = Inf)
  logLik <- -Inf
  while (any(err > unlist(options$tolerance[names(err)]))) {
    for (paramSet in names(fitSets)) {
      parNames <- fitSets[[paramSet]]
      res <- funs[[paramSet]](par)
      err[[paramSet]] <- getMaxAbsDifference(res, par[parNames])
      par[parNames] <- res
    }
    par["size"] <- fitParams(pulseData, par, "size", options)
    newlogLik <- evaluateLikelihood(par, pulseData)
    err["logLik"] <- - logLik + newlogLik
    logLik <- newlogLik
    log2screen(options, progressString(err, logLik))
    if (!is.null(options$resultRDS)) {
      ## assign names
      for(param in fitSets$params) {
        names(par[[param]]) <- rownames(pulseData$counts)
      }
      saveRDS(object = par, file = options$resultRDS)
    }
  }
  ## fit gene specific final parameters
  par[fitSets$params] <- funs[["params"]](par)
  if (!is.null(options$resultRDS)) {
    ## assign names
    for (param in fitSets$params) {
      names(par[[param]]) <- rownames(pulseData$counts)
    }
    saveRDS(object = par, file = options$resultRDS)
  }
  par
}

progressString <- function(err, logLik) {
  str <- format(unlist(err),
                digits = 2,
                width = 6)
  logLik <- format(logLik, digits=2, width=6)
  paste0("LogLik: [", logLik, " ] Max Rel.err. in ",
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
